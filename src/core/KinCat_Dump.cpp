#include "KinCat_Dump.hpp"
#include "KinCat_Util.hpp"

namespace KinCat {

template <typename DT>
Dump<DT>::Dump(const std::string &filename_dump, const std::string &filename_site, const Lattice<DT> &lattice)
    : _filename_dump(filename_dump), _filename_site(filename_site), _lattice(lattice), _ofs(), _lattice_host(),
      _is_first(false) {}

template <typename DT>
ordinal_type Dump<DT>::initialize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Dump::initialize", __FILE__, __LINE__, verbose);
  /// open output stream
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

  //_ofs.close();
  _is_first = true;
  snapshot(t, verbose);
  //_is_first = false;

  return 0;
}

template <typename DT>
ordinal_type Dump<DT>::finalize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Dump::finalize", __FILE__, __LINE__, verbose);

  const std::string indent("    ");

  /// close the file
  {
    _ofs << "\n" << indent << "]\n";
    _ofs << "}\n";

    _ofs.close();
  }

  /// last snapshot backup
  {
    _ofs.open(_filename_site, std::ios::out | std::ios::trunc);
    _is_first = true;
    _ofs << "{\n";
    _ofs << indent << "\"sites\": [\n";
    snapshot(t, verbose);
    _ofs << indent << "]\n";
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
                                const value_type_1d_view<site_type, host_device_type> sites,
                                const ordinal_type verbose) {
  const FunctionScope func_scope("Dump::snapshot", __FILE__, __LINE__, verbose);

  if (sites.extent(0) > 0) {
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

  return 0;
}

template <typename DT>
ordinal_type Dump<DT>::snapshot(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Dump::snapshot::batch", __FILE__, __LINE__, verbose);

  std::stringstream ss;
  ss << "Error: fails to open a dump file: " << _filename_dump;
  KINCAT_CHECK_ERROR(!_ofs.is_open(), ss.str().c_str());

  Kokkos::deep_copy(_lattice_host._sites, _lattice._sites);

  for (ordinal_type i = 0, iend = _lattice_host._sites.extent(0); i < iend; ++i) {
    const auto sites = Kokkos::subview(_lattice_host._sites, i, Kokkos::ALL());
    snapshot(i, t(i), sites, verbose);
  }

  return 0;
}

template <typename DT>
std::ostream &Dump<DT>::showMe(std::ostream &os, const std::string &label, const ordinal_type verbose) const {
  const std::string indent("  ");

  os << "-- Dump : " << label << "\n";
  os << indent << "-- filename dump : \"" << _filename_dump << "\"\n";
  os << indent << "-- filename site: \"" << _filename_site << "\"\n";

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
