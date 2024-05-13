#include "KinCat_Dump.hpp"
#include "KinCat_Util.hpp"
#include "KinCat_HDF5.hpp"

namespace KinCat {

template <typename DT>
Dump<DT>::Dump(const std::string &filename_dump, const std::string &filename_site, const Lattice<DT> &lattice, const bool useHDF5)
    : _filename_dump(filename_dump), _filename_site(filename_site), _lattice(lattice), _ofs(), _lattice_host(),
      _is_first(false), _useHDF5(useHDF5)
{}

template <typename DT>
ordinal_type Dump<DT>::initialize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Dump::initialize", __FILE__, __LINE__, verbose);

  //Check filename for HDF5 usage
  std::string hdf5_ext = ".hdf5";
  std::string json_ext = ".json";
  bool hdf5_file = _filename_dump.find(hdf5_ext) != std::string::npos;
  bool json_file = _filename_dump.find(json_ext) != std::string::npos;
  KINCAT_CHECK_ERROR(hdf5_file == json_file, "Dump output file extension not recognized.")
  if (!_useHDF5) { //Can't use HDF5
    KINCAT_CHECK_ERROR((hdf5_file), "Dump output file indicates HDF5 use. HDF5 is not enabled in this build.")
  } else { //Can use HDF5
    if (json_file) {
      _useHDF5 = false; //Use json output, even though HDF5 available. 
    }
  }

  /// create host mirror object of lattice
  _lattice.createMirrorObject(_lattice_host);

  /// create initial snapshot
  const ordinal_type n_species = _lattice_host._n_species;
  const ordinal_type n_samples = _lattice_host._sites.extent(0);
  const ordinal_type n_sites = _lattice_host._sites.extent(1);
  KINCAT_CHECK_ERROR(n_sites <= 0, "Error: sites are not created in the lattice");

  if (_useHDF5)
  {
#ifdef HAVE_HDF5
    std::cout << "creating hdf5 file : " << _filename_dump << '\n';
    _hdf5.Create(_filename_dump);
    _hdf5.Write("Dump", "number of species", n_species);
    _hdf5.Write("Dump", "coordinates.shape.dim0", n_sites);
    _hdf5.Write("Dump", "coordinates.shape.dim1", 2);
    
    ordinal_type s(0), k0, k1, ib;
    real_type coord_x, coord_y;

    std::vector<double> coordinates(n_sites * 2);
    
    for (s = 0; s < n_sites; ++s) {
      _lattice_host.getLatticeIndex(s, k0, k1, ib);
      _lattice_host.getCoordinates(k0, k1, ib, coordinates[s*2], coordinates[s*2+1]);
    }
    hid_t type = H5T_NATIVE_DOUBLE;
    _hdf5.Write("Dump", "coordinates.data", type, coordinates.size(), &coordinates[0]);
#else
#endif
  }
  else // _useHDF5 == false
  {
    /// open output stream
    std::cout << "creating json file : " << _filename_dump << '\n';
    _ofs.open(_filename_dump, std::ios::out | std::ios::trunc);

    std::stringstream ss;
    ss << "Error: fails to open a file: " << _filename_dump;
    KINCAT_CHECK_ERROR(!_ofs.is_open(), ss.str().c_str());

    /// create host mirror object of lattice
    _lattice.createMirrorObject(_lattice_host);

    /// create initial snapshot
    const ordinal_type n_species = _lattice_host._n_species;
    const ordinal_type n_samples = _lattice_host._sites.extent(0);
    const ordinal_type n_sites = _lattice_host._sites.extent(1);
    KINCAT_CHECK_ERROR(n_sites <= 0, "Error: sites are not created in the lattice");

    const std::string indent("    "), indent2(indent + indent);
    _ofs << "{\n";
    _ofs << indent << "\"samples\" : " << n_samples << ",\n";
    _ofs << indent << "\"number of species\": " << n_species << ",\n";
    _ofs << indent << "\"coordinates\": {\n";
    _ofs << indent2 << "\"shape\": [ " << n_sites << ", 2 ], \n";

    ordinal_type s(0), k0, k1, ib;
    real_type coord_x, coord_y;

    _lattice_host.getLatticeIndex(s, k0, k1, ib);
    _lattice_host.getCoordinates(k0, k1, ib, coord_x, coord_y);

    _ofs << indent2 << "\"data\": [ " << coord_x << ", " << coord_y;
    for (s = 1; s < n_sites; ++s) {
      _lattice_host.getLatticeIndex(s, k0, k1, ib);
      _lattice_host.getCoordinates(k0, k1, ib, coord_x, coord_y);
      _ofs << ", " << coord_x << ", " << coord_y;
    }
    _ofs << "]\n";
    _ofs << indent << "},\n";

    _ofs << indent << "\"sites\": [\n";

    // we leave _ofs open until finalize() is called
  }
  _is_first = true;
  snapshot(t, verbose);
  //_is_first = false;

  return 0;
}

template <typename DT>
ordinal_type Dump<DT>::finalize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Dump::finalize", __FILE__, __LINE__, verbose);

  const std::string indent("    ");

  if (_useHDF5)
  {
#ifdef HAVE_HDF5
    /// record the number of sites
    _hdf5.Write("Dump", "number of snapshots", _nextSnapshotID);
    
    /// close the file
    //_hdf5.Close(); //Why is this throwing an error? 
    
    /// last snapshot backup
    //{
    //  _hdf5.Create(_filename_site);
    //  _is_first = true;
    //  _nextSnapshotID = 0;
    //  snapshot(t, verbose);
    //  _hdf5.Write("Dump", "number of snapshots", _nextSnapshotID);
    //  _hdf5.Close();
    //}
#else
#endif
  }
  else //Not using HDF5
  {
    /// close the file
    {
      _ofs << "\n" << indent << "]\n";
      _ofs << "}\n";

      _ofs.close();
    }
  }
  /// last snapshot backup, create json, even if using HDF5. No HDF5 read-in yet. 
  {
    _ofs.open(_filename_site, std::ios::out | std::ios::trunc);
    std::stringstream ss;
    ss << "Error: fails to open a file: " << _filename_site;
    KINCAT_CHECK_ERROR(!_ofs.is_open(), ss.str().c_str());

    _is_first = true;
    _ofs << "{\n";
    _ofs << indent << "\"sites\": {\n";
    sitessnapshot(t, verbose);
    _ofs << indent << "}\n";
    _ofs << "\n}\n";
    //_is_first = false;
    _ofs.close();
  }

  /// delete host mirror object
  _lattice_host = Lattice<host_device_type>();

  return 0;
}

template <typename DT>
ordinal_type Dump<DT>::snapshot(const ordinal_type sid, const real_type t,
                                const value_type_1d_view<site_type, host_device_type> sites, const bool restart_flag,
                                const ordinal_type verbose) {
  const FunctionScope func_scope("Dump::snapshot", __FILE__, __LINE__, verbose);

  if (sites.extent(0) > 0) {
    if (_useHDF5)
    {
#ifdef HAVE_HDF5
      std::string DataSetName = std::string("snapshot.") + std::to_string(_nextSnapshotID);
      
      _hdf5.Write("Dump", DataSetName + ".sample", sid);
      _hdf5.Write("Dump", DataSetName + ".time",   t);
      _hdf5.WriteShortView1D("Dump", DataSetName + ".data", sites);
      
      _nextSnapshotID++;
#else
#endif
    }
    else // _useHDF5 == false
    {
      const std::string indent("    "), indent2(indent + indent), indent3(indent + indent2);
      if (!_is_first)
        _ofs << ", \n";
      _ofs << indent2 << "{ \n";
      _ofs << indent3 << "\"sample\": " << sid << ",\n";
      _ofs << indent3 << "\"time\": " << t << ",\n";
      _ofs << indent3 << "\"data\": [ " << ordinal_type(sites(0));
      for (ordinal_type i = 1, iend = sites.extent(0); i < iend; ++i) {
        _ofs << ", " << ordinal_type(sites(i));
      }
      _ofs << " ]\n";
      _ofs << indent2 << "}";
      _is_first = false;
    }
    if (_useHDF5 && restart_flag) // _useHDF5 == false
    {
      const std::string indent("    "), indent2(indent + indent), indent3(indent + indent2);
      if (!_is_first)
        _ofs << ", \n";
      _ofs << indent2 << "{ \n";
      _ofs << indent3 << "\"sample\": " << sid << ",\n";
      _ofs << indent3 << "\"time\": " << t << ",\n";
      _ofs << indent3 << "\"data\": [ " << ordinal_type(sites(0));
      for (ordinal_type i = 1, iend = sites.extent(0); i < iend; ++i) {
        _ofs << ", " << ordinal_type(sites(i));
      }
      _ofs << " ]\n";
      _ofs << indent2 << "}";
      _is_first = false;
    }
  }

  return 0;
}

template <typename DT>
ordinal_type Dump<DT>::snapshot(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Dump::snapshot::batch", __FILE__, __LINE__, verbose);

  std::stringstream ss;
  ss << "Error: fails to open a dump file: " << _filename_dump;
  if (_useHDF5)
  {
#ifdef HAVE_HDF5
    KINCAT_CHECK_ERROR(!_hdf5.IsOpen(), ss.str().c_str());
#else
#endif
  }
  else
  {
    KINCAT_CHECK_ERROR(!_ofs.is_open(), ss.str().c_str());
  }

  Kokkos::deep_copy(_lattice_host._sites, _lattice._sites);

  for (ordinal_type i = 0, iend = _lattice_host._sites.extent(0); i < iend; ++i) {
    const auto sites = Kokkos::subview(_lattice_host._sites, i, Kokkos::ALL());
    snapshot(i, t(i), sites, false, verbose);
  }

  return 0;
}

template <typename DT>
ordinal_type Dump<DT>::sitessnapshot(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Dump::sitessnapshot::batch", __FILE__, __LINE__, verbose);

  std::stringstream ss;
  ss << "Error: fails to open a dump file: " << _filename_dump;
  KINCAT_CHECK_ERROR(!_ofs.is_open(), ss.str().c_str());

  Kokkos::deep_copy(_lattice_host._sites, _lattice._sites);
  const auto sites = Kokkos::subview(_lattice_host._sites, 0, Kokkos::ALL());
  _ofs << "        \"0\": \n";
  snapshot(0, t(0), sites, true, verbose);

  for (ordinal_type i = 1, iend = _lattice_host._sites.extent(0); i < iend; ++i) {
    const auto sites = Kokkos::subview(_lattice_host._sites, i, Kokkos::ALL());
    _ofs << ",\n        \"" << std::to_string(i) << "\": ";
    _is_first = true;
    snapshot(i, t(i), sites, true, verbose);
  }
  _ofs << '\n';
  return 0;
}

template <typename DT>
std::ostream &Dump<DT>::showMe(std::ostream &os, const std::string &label, const ordinal_type verbose) const {
  const std::string indent("  ");

  os << "-- Dump : " << label << "\n";
  os << indent << "-- filename dump : \"" << _filename_dump << "\"\n";
  os << indent << "-- filename restart: \"" << _filename_site << "\"\n";

  return os;
}

/// eti for individual device type
#if defined(KOKKOS_ENABLE_SERIAL)
template struct Dump<typename UseThisDevice<Kokkos::Serial>::type>;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template struct Dump<typename UseThisDevice<Kokkos::OpenMP>::type>;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template struct Dump<typename UseThisDevice<Kokkos::Cuda>::type>;
#endif

} // namespace KinCat
