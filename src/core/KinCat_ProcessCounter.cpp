#include "KinCat_ProcessCounter.hpp"

namespace KinCat {

template <typename DT>

ProcessCounter<DT>::ProcessCounter(const ordinal_type n_domains, const ordinal_type n_processes,
                                   const value_type_2d_view<ordinal_type, DT> &processints)
    : _is_initialized(false), _processints(processints),_n_domains(n_domains), _n_processes(n_processes), _process_counter(),
      _processints_counter(), _process_counter_host(), _processints_counter_host() {}


template <typename DT> void ProcessCounter<DT>::initialize(const ordinal_type n_samples) {
  _process_counter = value_type_2d_view<size_type, DT>("process counter", n_samples, _n_processes*_n_domains);
  _n_processints = _processints.extent(0);
  _processints_counter =
      value_type_2d_view<size_type, DT>("process instance counter", n_samples, (_processints.extent(0)*_n_domains) + 1);
  _process_counter_host = Kokkos::create_mirror_view(host_device_type(), _process_counter);
  _processints_counter_host = Kokkos::create_mirror_view(host_device_type(), _processints_counter);

  //_functiontimer = value_type_1d_view<std::chrono::duration<double>, DT>("function timers", 5);

  _is_initialized = true;
}

template <typename DT> void ProcessCounter<DT>::reset() const {
  const size_type zero(0);
  Kokkos::deep_copy(_process_counter, zero);
  Kokkos::deep_copy(_processints_counter, zero);
}

template <typename DT> void ProcessCounter<DT>::syncToHost() const {
  Kokkos::deep_copy(_process_counter_host, _process_counter);
  Kokkos::deep_copy(_processints_counter_host, _processints_counter);
}

template <typename DT>
value_type_1d_view<size_type, host_device_type>
ProcessCounter<DT>::getProcessCounterHost(const ordinal_type sid) const {
  return Kokkos::subview(_process_counter_host, sid, Kokkos::ALL());
}

template <typename DT>
value_type_1d_view<size_type, host_device_type>
ProcessCounter<DT>::getProcessIntsCounterHost(const ordinal_type sid) const {
  return Kokkos::subview(_processints_counter_host, sid, Kokkos::ALL());
}

template <typename DT>
ordinal_type
ProcessCounter<DT>::getTotalProcessCount(const ordinal_type s) const {
  return _processints_counter_host(s, _processints_counter_host.extent(1) - 1);
}

template <typename DT>
std::ostream &ProcessCounter<DT>::showMe(std::ostream &os, const std::string &label, const ordinal_type verbose) const {
  const std::string indent2("  "), indent4("  "), indet6("      ");
  //os << "n_domains in Counter : " << _n_domains << '\n';
  // os << "-- ProcessCounter : " << label << "\n";
  // os << indent2 << "-- is initialized : " << (_is_initialized ? "true" : "false") << "\n";

  if (_is_initialized) {
    syncToHost();
    const ordinal_type n_samples = _process_counter_host.extent(0);
    for (ordinal_type s = 0, send = n_samples; s < send; ++s) {
      os << indent2
         << "-- # of events occured : " << _processints_counter_host(s, _n_domains*_n_processints) << "\n"; //(s, _processints_counter_host.extent(1) - 1) ???
      os << indent2 << "-- sample id : " << s << "\n";
      os << indent2 << "-- process counter : \n";
      for (ordinal_type i = 0, iend = _process_counter_host.extent(1); i < iend; ++i) {
        ordinal_type tmp_counter = 0;
        for (ordinal_type k = 0, kend = _n_domains; k < kend; ++k) {
          tmp_counter += _processints_counter_host(s, i + (k * iend));
        }
        os << indent4 << "[ " << i << ", " << tmp_counter << "] \n";
      }

      if (verbose > 3) {
        os << indent2 << "-- process instance counter : \n";
        ordinal_type tmp_counter;
        for (ordinal_type i = 0, iend = (_n_processints); i < iend; ++i) {
          tmp_counter = 0;
          for (ordinal_type k = 0, kend = _n_domains; k < kend; ++k) {
            tmp_counter += _processints_counter_host(s, i + (k*iend));
          }
          os << indent4 << "[ " << i << ", " << tmp_counter << "] \n";
        }
      }
    }
  }
  return os;
}

template <typename DT>
std::ostream &ProcessCounter<DT>::showMe(std::ostream &os, const std::string &label,
                                         const std::vector<std::string> &processes, const ordinal_type verbose) const {
  const std::string indent2("  "), indent4("  "), indet6("      ");
  if (_is_initialized) {
    syncToHost();
    const ordinal_type n_samples = _process_counter_host.extent(0);
    for (ordinal_type s = 0, send = n_samples; s < send; ++s) {
      if (n_samples > 1) {
        os << indent2 << "-- sample id : " << s << "\n";
      }
      os << indent2
         //<< "-- # of events occured : " << _processints_counter_host(s, _processints_counter_host.extent(1) - 1)
         << "-- # of events occured : " << _processints_counter_host(s, _n_domains*_n_processints) //1 ???
         << "\n";
      //std::cout << _processints_counter_host.extent(0) << ", " << _processints_counter_host.extent(1) << '\n';
      if (verbose > 3) {
        os << indent2 << "-- process instance counter : \n";
        ordinal_type tmp_counter;
        for (ordinal_type i = 0, iend = (_n_processints); i < iend; ++i) {
          tmp_counter = 0;
          for (ordinal_type j = 0, jend = _n_domains; j < jend; ++j) {
            tmp_counter += _processints_counter_host(s, i+(j*iend));
          }
          os << indent4 << "[ " << i << ", " << tmp_counter << "] \n";
        } 
      }
    }
  }
  return os;
}

/// eti for individual device type
#if defined(KOKKOS_ENABLE_SERIAL)
template struct ProcessCounter<typename UseThisDevice<Kokkos::Serial>::type>;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template struct ProcessCounter<typename UseThisDevice<Kokkos::OpenMP>::type>;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template struct ProcessCounter<typename UseThisDevice<Kokkos::Cuda>::type>;
#endif

} // namespace KinCat
