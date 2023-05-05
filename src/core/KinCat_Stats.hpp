#ifndef __KINCAT_STATS_HPP__
#define __KINCAT_STATS_HPP__

#include "KinCat_Lattice.hpp"
#include "KinCat_ProcessCounter.hpp"
#include "KinCat_Util.hpp"

namespace KinCat {

template <typename DeviceType> struct Stats {
private:
  std::string _filename_stats;
  std::vector<std::string> _stats_list;
  std::vector<std::string> _processes;
  ordinal_type _n_processes;
  ordinal_type _n_basis_sites; 
  Lattice<DeviceType> _lattice;

  std::ofstream _ofstats;
  Lattice<host_device_type> _lattice_host;
  ProcessCounter<DeviceType> _counter;

  // Statistic type flags
  bool _s_c; // species coverage
  bool _p_c; // process counts
  bool _s_s_c; // site specific species coverage

  bool _is_first;

  ordinal_type snapshot(const ordinal_type sid, const real_type t,
                        const value_type_1d_view<site_type, host_device_type> sites, const ordinal_type verbose = 0);

public:
  Stats() = delete;
  Stats(const Stats &b) = delete;
  Stats(const std::string &filename_stats, const std::vector<std::string> &stats_list,
        const Lattice<DeviceType> &lattice, const std::vector<std::string> &processes,
        const ProcessCounter<DeviceType> &counter);
  virtual ~Stats() = default;

  ordinal_type initialize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);
  ordinal_type finalize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);
  ordinal_type snapshot(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);

  virtual std::ostream &showMe(std::ostream &os, const std::string &label, const ordinal_type verbose = 0) const;
};

} // namespace KinCat

#endif
