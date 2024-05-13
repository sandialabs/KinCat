#include "KinCat_Stats.hpp"
#include "KinCat_Util.hpp"
#include "KinCat_HDF5.hpp"

namespace KinCat {

template <typename DT>
Stats<DT>::Stats(const std::string &filename_stats, const std::vector<std::string> &stats_list,
                 const Lattice<DT> &lattice, const std::vector<std::string> &processes,
                 const ProcessCounter<DT> &counter, const bool useHDF5)
    : _filename_stats(filename_stats), _stats_list(stats_list), _lattice(lattice), _ofstats(), _lattice_host(),
      _is_first(false), _processes(processes), _n_processes(processes.size()), _useHDF5(useHDF5) {
  _counter = counter;
  _n_domains = counter.getProcessCounterHost().extent(0)/_n_processes;
}

template <typename DT>
ordinal_type Stats<DT>::initialize(const value_type_1d_view<real_type, host_device_type> t,
                                   const ordinal_type verbose) {
  const FunctionScope func_scope("Stats::initialize", __FILE__, __LINE__, verbose);
  
  //Check filename for HDF5 usage
  std::string hdf5_ext = ".hdf5";
  std::string json_ext = ".json";
  bool hdf5_file = _filename_stats.find(hdf5_ext) != std::string::npos;
  bool json_file = _filename_stats.find(json_ext) != std::string::npos;
  KINCAT_CHECK_ERROR(hdf5_file == json_file, "Stats output file extension not recognized.")
  if (!_useHDF5) { //Can't use HDF5
    KINCAT_CHECK_ERROR((hdf5_file), "Stats output file indicates HDF5 use. HDF5 is not enabled in this build.")
  } else { //Can use HDF5
    if (json_file) {
      _useHDF5 = false; //Use json output, even though HDF5 available. 
    }
  }
  

  /// create host mirror object of lattice
  _lattice.createMirrorObject(_lattice_host);

  /// determine which statistics are to be included, Add new statistics names here.
  std::unordered_map<std::string, std::string> stat_dictionary; // set of possible/recognized keys
  stat_dictionary["species_coverage"] = "species_coverage";
  stat_dictionary["process_counts"] = "process_counts";
  stat_dictionary["site_species_coverage"] = "site_species_coverage";

  std::unordered_map<std::string, std::string> stat_keys; // set of keys included by user
  for (ordinal_type i = 0, iend = _stats_list.size(); i < iend; ++i) {
    stat_keys[_stats_list[i]] = _stats_list[i];
    bool key_flag = (!stat_dictionary.count(_stats_list[i]));
    std::stringstream ee;
    ee << "Unrecognized type in statistics: " << _stats_list[i];
    KINCAT_CHECK_ERROR(key_flag, ee.str().c_str());
  }
  // Initialize flags for statistics types (flags should be instantiated in KinCat_Stats.hpp)
  if (stat_keys.count("species_coverage")) {
    _s_c = true;
  }
  if (stat_keys.count("process_counts")) {
    _p_c = true;
  }
  if (stat_keys.count("site_species_coverage")) {
    _s_s_c = true;
  }

  /// create initial snapshot
  const ordinal_type n_species = _lattice_host._n_species;
  const ordinal_type n_samples = _lattice_host._sites.extent(0);
  const ordinal_type n_sites = _lattice_host._sites.extent(1);
  const ordinal_type n_cells_x = _lattice_host._n_cells_x;
  const ordinal_type n_cells_y = _lattice_host._n_cells_y;
  const ordinal_type n_basis_sites = n_sites / (n_cells_x * n_cells_y);
  _n_basis_sites = n_basis_sites; //set class variable so can access in snapshot
  KINCAT_CHECK_ERROR(n_sites <= 0, "Error: sites are not created in the lattice");

  if (_useHDF5) {
#ifdef HAVE_HDF5
    std::cout << "creating hdf5 stats file : " << _filename_stats << '\n';
    _hdf5_stats.Create(_filename_stats);
    _hdf5_stats.Write("Stats", "lattice.shape.dim1", n_cells_x);
    _hdf5_stats.Write("Stats", "lattice.shape.dim2", n_cells_y);
    _hdf5_stats.Write("Stats", "lattice.shape.basis", n_basis_sites);
    //_hdf5_stats.Write("Stats", "samples", n_samples);
    //_hdf5_stats.Write("Stats", "number.of.species", n_species);
    //_hdf5_stats.Write("Stats", "number.of.processes", _n_processes);
    for (ordinal_type i = 0, iend = _n_processes; i < iend; i++) {
      std::string ProcessName = std::string("process.") + std::to_string(i);
      //_hdf5_stats.Write("Stats", ProcessName, _processes[i]);
    }

#else
#endif
  }
  else { // _useHDF5 == false
    /// open output stream
    _ofstats.open(_filename_stats, std::ios::out | std::ios::trunc);

    std::stringstream ss;
    ss << "Error: fails to open a statistics file: " << _filename_stats;
    KINCAT_CHECK_ERROR(!_ofstats.is_open(), ss.str().c_str());

    const std::string indent("    "), indent2(indent + indent);
    _ofstats << "{\n";
    // Initialize file header/formating of data.
    _ofstats << indent << "\"lattice size\" : [ " << n_cells_x << ", " << n_cells_y << ", " << n_basis_sites << " ], \n";
    _ofstats << indent << "\"samples\" : " << n_samples << ",\n";
    if (_s_c || _s_s_c) {
      _ofstats << indent << "\"number of species\" : " << n_species << ",\n";
    }
    if (_p_c) {
      _ofstats << indent << "\"number of processes\" : " << _n_processes << ",\n";
      _ofstats << indent << "\"processes\" : [ \"" << _processes[0] << "\"";
      for (ordinal_type i = 1; i < _n_processes; i++) {
        _ofstats << ", \"" << _processes[i] << "\"";
      }
      _ofstats << " ], \n";
    }

    _ofstats << indent << "\"readings\" : [ \n";

    _is_first = true;
    snapshot(t, verbose);
    //_is_first = false; // now done in snapshot
  }

  return 0;
}

template <typename DT>
ordinal_type Stats<DT>::finalize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Stats::finalize", __FILE__, __LINE__, verbose);

  const std::string indent("    ");

  if (_useHDF5) {
#ifdef HAVE_HDF5
    _hdf5_stats.Write("Stats", "number of snapshots", _nextSnapshotID);
#else
#endif
  }
  else //Not using HDF5
  {  /// close the file
    {
      _ofstats << "\n" << indent << "]\n";
      _ofstats << "}\n";

      _ofstats.close();
    }
  }
  
  /// delete host mirror object
  _lattice_host = Lattice<host_device_type>();

  return 0;
}

template <typename DT>
ordinal_type Stats<DT>::snapshot(const ordinal_type sid, const real_type t,
                                 const value_type_1d_view<site_type, host_device_type> sites,
                                 const ordinal_type verbose) {
  const FunctionScope func_scope("Stats::snapshot", __FILE__, __LINE__, verbose);

  std::stringstream ss;
  ss << "Error: fails to open a file: " << _filename_stats;
  const ordinal_type n_species = _lattice_host._n_species;
  const ordinal_type n_sites = sites.extent(0);

  if (_useHDF5) {
#ifdef HAVE_HDF5

    KINCAT_CHECK_ERROR(!_hdf5_stats.IsOpen(), ss.str().c_str());
    if (sites.extent(0) > 0) {
      std::string DataSetName = std::string("snapshot.") + std::to_string(_nextSnapshotID);
      _hdf5_stats.Write("Stats", DataSetName + ".sample", sid);
      _hdf5_stats.Write("Stats", DataSetName + ".time", t);

      if (_s_c) {
        //_ofstats << ",\n" << indent3 << "\"species coverage\" : [ ";
        ordinal_type spec_coverages[n_species];
        for (ordinal_type i = 0, iend = n_species; i < iend; ++i) {
          spec_coverages[i] = 0;
        }
        for (ordinal_type i = 0, iend = n_sites; i < iend; ++i) {
          spec_coverages[ordinal_type(sites(i))]++;
        }
        real_type denom; //coverage denominator
        denom = sites.size();

        value_type_1d_view<real_type, host_device_type> coverage ("coverage view", n_species);

        //real_type coverage = real_type(spec_coverages[0]) / denom;
        //_ofstats << coverage;
        for (ordinal_type i = 0, iend = n_species; i < iend; ++i) {
          coverage(i) = real_type(spec_coverages[i]) / denom;
          //_ofstats << ", " << coverage;
        }
        _hdf5_stats.WriteDoubleView1D("Stats", DataSetName + ".coverage", coverage);
        //_ofstats << " ]";
      }
      if (_p_c) {
        value_type_1d_view<size_type, host_device_type> process_count = _counter.getProcessCounterHost(sid);
        value_type_1d_view<int, host_device_type> process_count_ints ("process counter ints", process_count.extent(0));
        ordinal_type temp_count = 0;
        for (ordinal_type i = 0, iend = process_count.extent(0); i < iend; i++) {
          temp_count = 0;
          for (ordinal_type j = 0, jend = _n_domains; j < jend; j++) {
            temp_count += int(process_count(i + (j*_n_processes)));
          }
          process_count_ints(i) = temp_count;//int(process_count(i));
        }
        _hdf5_stats.WriteIntView1D("Stats", DataSetName + ".processcounts", process_count_ints);
      }
      if (_s_s_c) {
        //_ofstats << ",\n" << indent3 << "\"site species coverage\" : [[ ";
        ordinal_type spec_coverages[_n_basis_sites][n_species];
        for (ordinal_type j = 0, jend = _n_basis_sites; j < jend; j++) {
          for (ordinal_type i = 0, iend = n_species; i < iend; ++i) {
            spec_coverages[j][i] = 0;
          }
        }
        for (ordinal_type i = 0, iend = n_sites; i < iend; ++i) {
          ordinal_type site_type = i % _n_basis_sites;
          spec_coverages[site_type][ordinal_type(sites(i))]++;
        }
        real_type denom; //coverage denominator
        value_type_1d_view<double, host_device_type> spec_coverages_view ("species coverages view", (_n_basis_sites*n_species));
        denom = sites.size()/_n_basis_sites; //equivalent to n_cells_x * n_cells_y
        real_type coverage;
        for (ordinal_type i = 0, iend = _n_basis_sites; i < iend; ++i) {
          for (ordinal_type j = 0, jend = n_species; j < jend; ++j) {
            coverage = real_type(spec_coverages[i][j]) / denom;
            spec_coverages_view((i*jend)+j) = coverage;
          }
        }
        _hdf5_stats.WriteDoubleView1D("Stats", DataSetName + ".speciescoverage", spec_coverages_view);
      }
      _is_first = false;

      _nextSnapshotID++;

    }
#else
#endif
  } else {
    KINCAT_CHECK_ERROR(!_ofstats.is_open(), ss.str().c_str());
    if (sites.extent(0) > 0) {
      const std::string indent("    "), indent2(indent + indent), indent3(indent + indent2);
      if (!_is_first) {
        _ofstats << ", \n";
      }
      _ofstats << indent2 << "{ \n";
      _ofstats << indent3 << "\"sample\" : " << sid << ",\n";
      _ofstats << indent3 << "\"time\" : " << t;
      // Output statistics
      if (_s_c) {
        _ofstats << ",\n" << indent3 << "\"species coverage\" : [ ";
        ordinal_type spec_coverages[n_species];
        for (ordinal_type i = 0, iend = n_species; i < iend; ++i) {
          spec_coverages[i] = 0;
        }
        for (ordinal_type i = 0, iend = n_sites; i < iend; ++i) {
          spec_coverages[ordinal_type(sites(i))]++;
        }
        real_type denom; //coverage denominator
        denom = sites.size();
        real_type coverage = real_type(spec_coverages[0]) / denom;
        _ofstats << coverage;
        for (ordinal_type i = 1, iend = n_species; i < iend; ++i) {
          coverage = real_type(spec_coverages[i]) / denom;
          _ofstats << ", " << coverage;
        }
        _ofstats << " ]";
      }
      if (_p_c) {
        value_type_1d_view<size_type, host_device_type> process_count = _counter.getProcessCounterHost(sid);
        
        ordinal_type temp_count = 0;
        for (ordinal_type j = 0, jend = _n_domains; j < jend; ++j) {
          temp_count += int(process_count(j*_n_processes));
        }
        _ofstats << ",\n" << indent3 << "\"process counts\" : [ " << temp_count;
        for (ordinal_type i = 1, iend = _n_processes; i < iend; ++i) {
          temp_count = 0;
          for (ordinal_type j = 0, jend = _n_domains; j < jend; ++j) {
            temp_count += int(process_count(i + j*_n_processes));
          }
          _ofstats << ", " << temp_count;
        }
        _ofstats << " ]";
      }
      if (_s_s_c) {
        _ofstats << ",\n" << indent3 << "\"site species coverage\" : [[ ";
        ordinal_type spec_coverages[_n_basis_sites][n_species];
        for (ordinal_type j = 0, jend = _n_basis_sites; j < jend; j++) {
          for (ordinal_type i = 0, iend = n_species; i < iend; ++i) {
            spec_coverages[j][i] = 0;
          }
        }
        for (ordinal_type i = 0, iend = n_sites; i < iend; ++i) {
          ordinal_type site_type = i % _n_basis_sites;
          spec_coverages[site_type][ordinal_type(sites(i))]++;
        }
        real_type denom; //coverage denominator
        denom = sites.size()/_n_basis_sites; //equivalent to n_cells_x * n_cells_y
        for (ordinal_type i = 0, iend = _n_basis_sites; i < iend; ++i) {
          real_type coverage = real_type(spec_coverages[i][0]) / denom;
          _ofstats << coverage;
          for (ordinal_type j = 1, jend = n_species; j < jend; ++j) {
            coverage = real_type(spec_coverages[i][j]) / denom;
            _ofstats << ", " << coverage;
          }
          if (i < iend-1) {
            _ofstats << " ], [ ";
          }
          else {
            _ofstats << " ]";
          }
        }
        _ofstats << "]";
      }
      _ofstats << '\n' << indent2 << "}";
      _is_first = false;
    }
  }

  return 0;
}

template <typename DT>
ordinal_type Stats<DT>::snapshot(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Stats::snapshot::batch", __FILE__, __LINE__, verbose);

  std::stringstream ss;
  ss << "Error: fails to open a file: " << _filename_stats;
  if (_useHDF5) {
#ifdef HAVE_HDF5
    KINCAT_CHECK_ERROR(!_hdf5_stats.IsOpen(), ss.str().c_str());
#else
#endif
  } else { 
    KINCAT_CHECK_ERROR(!_ofstats.is_open(), ss.str().c_str());
  }

  Kokkos::deep_copy(_lattice_host._sites, _lattice._sites); // Needed?

  for (ordinal_type i = 0, iend = _lattice_host._sites.extent(0); i < iend; ++i) {
    const auto sites = Kokkos::subview(_lattice_host._sites, i, Kokkos::ALL());
    snapshot(i, t(i), sites, verbose);
  }

  return 0;
}

template <typename DT>
std::ostream &Stats<DT>::showMe(std::ostream &os, const std::string &label, const ordinal_type verbose) const {
  const std::string indent("  ");

  os << "-- Stats : " << label << "\n";
  os << indent << "-- filename stats : \"" << _filename_stats << "\"\n";
  
  return os;
}

/// eti for individual device type
#if defined(KOKKOS_ENABLE_SERIAL)
template struct Stats<typename UseThisDevice<Kokkos::Serial>::type>;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template struct Stats<typename UseThisDevice<Kokkos::OpenMP>::type>;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template struct Stats<typename UseThisDevice<Kokkos::Cuda>::type>;
#endif

} // namespace KinCat
