#ifndef __KINCAT_PROCESS_COUNTER_HPP__
#define __KINCAT_PROCESS_COUNTER_HPP__

#include "KinCat_Util.hpp"

namespace KinCat {

template <typename DeviceType> struct ProcessCounter {
public:
  bool _is_initialized;

  value_type_2d_view<ordinal_type, DeviceType> _processints;
  ordinal_type _n_processes;
  ordinal_type _n_processints;
  ordinal_type _n_domains;

  /// counter
  value_type_2d_view<size_type, DeviceType> _process_counter;     /// (n_samples, n_processes * n_domains)
  value_type_2d_view<size_type, DeviceType> _processints_counter; /// (n_samples, (n_processints * n_domains) + 1)

  /// counter mirror to host
  value_type_2d_view<size_type, host_device_type> _process_counter_host;
  value_type_2d_view<size_type, host_device_type> _processints_counter_host;

  value_type_3d_view<size_type, host_device_type> _message_data_host;

  /// timers This functionality is currently removed, slows code. 
  //value_type_1d_view<std::chrono::duration<double>, DeviceType> _functiontimer; 

public:
  KOKKOS_DEFAULTED_FUNCTION ProcessCounter() = default;
  KOKKOS_DEFAULTED_FUNCTION ProcessCounter(const ProcessCounter &b) = default;
  KOKKOS_DEFAULTED_FUNCTION ~ProcessCounter() = default;
  ProcessCounter(const ordinal_type n_domains, const ordinal_type n_processes, const value_type_2d_view<ordinal_type, DeviceType> &processints);

  void initialize(const ordinal_type n_samples = 1);
  void reset() const;
  void syncToHost() const;

  value_type_1d_view<size_type, host_device_type> getProcessCounterHost(const ordinal_type sid = 0) const;
  value_type_1d_view<size_type, host_device_type> getProcessIntsCounterHost(const ordinal_type sid = 0) const;
  ordinal_type getTotalProcessCount(const ordinal_type s) const;

  KOKKOS_INLINE_FUNCTION
  void update(const ordinal_type sid, const ordinal_type did, const ordinal_type eid) const {
    Kokkos::atomic_add(&_process_counter(sid, _processints(eid, 2) + (did*_n_processes)), 1);
    Kokkos::atomic_add(&_processints_counter(sid, eid + (did*_n_processints)) , 1);

    /// optional count for processints occured
    const ordinal_type last = _processints_counter.extent(1) - 1;
    Kokkos::atomic_add(&_processints_counter(sid, last), 1);
  }

  //Timing functionality currently removed. 
  //void update_timer(const ordinal_type fid, const std::chrono::duration<double> t_incr) {
  //  _functiontimer(fid) += t_incr;
  //}

  virtual std::ostream &showMe(std::ostream &os, const std::string &label, const ordinal_type verbose = 0) const;
  std::ostream &showMe(std::ostream &os, const std::string &label, const std::vector<std::string> &processes,
                       const ordinal_type verbose = 0) const;
};

} // namespace KinCat

#endif
