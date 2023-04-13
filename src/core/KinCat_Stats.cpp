#include "KinCat_Stats.hpp"
#include "KinCat_Util.hpp"

namespace KinCat {

template <typename DT>
Stats<DT>::Stats(const std::string &filename_stats, const std::vector<std::string> &stats_list,
                 const Lattice<DT> &lattice, const std::vector<std::string> &processes,
                 const ProcessCounter<DT> &counter)
    : _filename_stats(filename_stats), _stats_list(stats_list), _lattice(lattice), _ofstats(), _lattice_host(),
      _is_first(false), _processes(processes), _n_processes(processes.size()) {
  _counter = counter;
}

template <typename DT>
ordinal_type Stats<DT>::initialize(const value_type_1d_view<real_type, host_device_type> t,
                                   const ordinal_type verbose) {
  const FunctionScope func_scope("Stats::initialize", __FILE__, __LINE__, verbose);
  /// open output stream
  _ofstats.open(_filename_stats, std::ios::out | std::ios::trunc);

  std::stringstream ss;
  ss << "Error: fails to open a statistics file: " << _filename_stats;
  KINCAT_CHECK_ERROR(!_ofstats.is_open(), ss.str().c_str());

  /// create host mirror object of lattice
  _lattice.createMirrorObject(_lattice_host);

  /// determine which statistics are to be included, Add new statistics names here.
  std::unordered_map<std::string, std::string> stat_dictionary; // set of possible/recognized keys
  stat_dictionary["species_coverage"] = "species_coverage";
  stat_dictionary["process_counts"] = "process_counts";

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
    // std::cout << "found species_coverage \n";
  }
  if (stat_keys.count("process_counts")) {
    _p_c = true;
  }

  /// create initial snapshot
  const ordinal_type n_species = _lattice_host._n_species;
  const ordinal_type n_samples = _lattice_host._sites.extent(0);
  const ordinal_type n_sites = _lattice_host._sites.extent(1);
  const ordinal_type n_cells_x = _lattice_host._n_cells_x;
  const ordinal_type n_cells_y = _lattice_host._n_cells_y;
  const ordinal_type n_basis_sites = n_sites / (n_cells_x * n_cells_y);
  KINCAT_CHECK_ERROR(n_sites <= 0, "Error: sites are not created in the lattice");

  const std::string indent("    "), indent2(indent + indent);
  _ofstats << "{\n";
  // Initialize file header/formating of data.
  _ofstats << indent << "\"lattice size\" : [ " << n_cells_x << ", " << n_cells_y << ", " << n_basis_sites << " ], \n";
  _ofstats << indent << "\"samples\" : " << n_samples << ",\n";
  if (_s_c) {
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

  return 0;
}

template <typename DT>
ordinal_type Stats<DT>::finalize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Stats::finalize", __FILE__, __LINE__, verbose);

  const std::string indent("    ");

  /// close the file
  {
    _ofstats << "\n" << indent << "]\n";
    _ofstats << "}\n";

    _ofstats.close();
  }
  //
  //  /// last snapshot backup
  //  {
  //    _ofstats.open(_filename_site, std::ios::out | std::ios::trunc);
  //    _is_first = true;
  //    _ofstats << "{\n";
  //    _ofstats << indent << "\"sites\": [\n";
  //    snapshot(t, verbose);
  //    _ofstats << indent << "]\n";
  //    _ofstats << "\n}\n";
  //    //_is_first = false;
  //    _ofstats.close();
  //  }
  //
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
  KINCAT_CHECK_ERROR(!_ofstats.is_open(), ss.str().c_str());

  // const auto sites = _lattice_host._sites;
  // const auto sites = Kokkos::subview(_lattice_host._sites, i, Kokkos::ALL());
  const ordinal_type n_species = _lattice_host._n_species;
  const ordinal_type n_sites = sites.extent(0);
  // Kokkos::deep_copy(sites, _lattice._sites);

  if (sites.extent(0) > 0) {
    const std::string indent("    "), indent2(indent + indent), indent3(indent + indent2);
    if (!_is_first)
      _ofstats << ", \n";
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
      real_type coverage = real_type(spec_coverages[0]) / sites.size();
      _ofstats << coverage;
      for (ordinal_type i = 1, iend = n_species; i < iend; ++i) {
        coverage = real_type(spec_coverages[i]) / sites.size();
        _ofstats << ", " << coverage;
      }
      _ofstats << " ]";
    }
    if (_p_c) {
      value_type_1d_view<size_type, host_device_type> process_count = _counter.getProcessCounterHost(sid);
      _ofstats << ",\n" << indent3 << "\"process counts\" : [ " << process_count(0);
      for (ordinal_type i = 1, iend = _n_processes; i < iend; ++i) {
        _ofstats << ", " << process_count(i);
      }
      _ofstats << " ]";
    }
    _ofstats << '\n' << indent2 << "}";
    _is_first = false;
  }

  return 0;
}

template <typename DT>
ordinal_type Stats<DT>::snapshot(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose) {
  const FunctionScope func_scope("Stats::snapshot::batch", __FILE__, __LINE__, verbose);

  std::stringstream ss;
  ss << "Error: fails to open a file: " << _filename_stats;
  KINCAT_CHECK_ERROR(!_ofstats.is_open(), ss.str().c_str());

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
  // os << indent << "-- filename site: \"" << _filename_site << "\"\n";

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
