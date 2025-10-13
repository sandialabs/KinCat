// clang-format off
#include "KinCat_Util.hpp"
#include "KinCat_SolverSerialRTA.hpp"
#include "KinCat_SolverInternal.hpp"
// clang-format on

namespace KinCat {

template <typename DT>
SolverSerialRTA<DT>::SolverSerialRTA() : _sid(0), _instance_rates(), _instance_rates_scan(), _random_pool() {
  const bool device_has_host_memory = std::is_same<typename DT::memory_space, Kokkos::HostSpace>::value;
  KINCAT_CHECK_ERROR(!device_has_host_memory,
                     "Error: SolverSerialRTA can only be used with host memory space i.e., no CUDA");
}
template <typename DT> SolverSerialRTA<DT>::SolverSerialRTA(const SolverSerialRTA &b) = default;
template <typename DT> SolverSerialRTA<DT>::~SolverSerialRTA() = default;

template <typename DT> ordinal_type SolverSerialRTA<DT>::fillRandomPool(const ordinal_type verbose) {
  const FunctionScope func_scope("SolverSerialRTA::fillRandomPool", __FILE__, __LINE__, verbose);
  Kokkos::fill_random(_random_pool, _random, real_type(1));
  _iter_random_pool = 0;
  if (verbose > 2) {
    const ordinal_type iend = std::min<ordinal_type>(_random_pool.extent(0), 10);
    std::cout << "random pool\n";
    std::cout << "  [ " << _random_pool(0);
    for (ordinal_type i = 1; i < iend; ++i) {
      std::cout << ", " << _random_pool(i);
    }
    std::cout << "]\n";
  }
  return 0;
}

template <typename DT> real_type SolverSerialRTA<DT>::getRandomNumber(const ordinal_type verbose) {
  const FunctionScope func_scope("SolverSerialRTA::getRandomNumber", __FILE__, __LINE__, verbose);
  if (_iter_random_pool >= ordinal_type(_random_pool.extent(0))) {
    fillRandomPool(verbose);
  }
  const real_type r_val = _random_pool[_iter_random_pool++];
  if (verbose > 2) {
    std::cout << "random value: " << r_val << "\n";
  }
  return r_val;
}

template <typename DT>
ordinal_type
SolverSerialRTA<DT>::initialize(const bool is_solver_random_number_variation, const ordinal_type random_seed,
                                const ordinal_type random_pool_size, const ordinal_type max_kmc_steps_per_kernel_launch,
                                const ordinal_type n_cells_interaction_x, const ordinal_type n_cells_interaction_y,
                                const Lattice<DT> &arg_lattice, const ProcessDictionary<DT> &arg_dictionary,
                                const ProcessCounter<DT> &arg_counter, const ordinal_type verbose) {
  const FunctionScope func_scope("SolverSerialRTA::initialize", __FILE__, __LINE__, verbose);

  _sid = 0;

  const bool device_has_host_memory = std::is_same<typename DT::memory_space, Kokkos::HostSpace>::value;
  KINCAT_CHECK_ERROR(!device_has_host_memory, "Error: serial-rta solver can only be used with host memory space i.e., "
                                              "no CUDA\n Use batch-rta with a single sample instead");

  /// set member variables
  _lattice = arg_lattice;
  _dictionary = arg_dictionary;
  _counter = arg_counter;

  /// serial version always has the same domain as the lattice
  _lattice.setDomain(_lattice._n_cells_x, _lattice._n_cells_y);

  /// set random variable
  _iter_random_pool = 0;
  _random = Kokkos::Random_XorShift64_Pool<DT>(random_seed);

  /// set max kmc steps per kernel launch (used for workspace and random pool setup)
  _max_kmc_steps_per_kernel_launch = max_kmc_steps_per_kernel_launch;

  /// uniform random values (we do not want to refill random numbers during a kmc run)
  const ordinal_type required_random_pool_size = _max_kmc_steps_per_kernel_launch * 2;
  const ordinal_type p0 = random_pool_size / required_random_pool_size;
  const ordinal_type p1 = random_pool_size % required_random_pool_size;
  const ordinal_type actual_random_pool_size = (p0 + (p1 > 0)) * required_random_pool_size;

  _random_pool = value_type_1d_view<real_type, DT>(do_not_init_tag("random pool"), actual_random_pool_size);
  KINCAT_CHECK_ERROR(fillRandomPool(verbose), "Error: fillRandomPool returns non-zero error code");

  /// create instance list
  const ordinal_type n_cells = _lattice.getNumberOfCells();
  _instance_rates = value_type_1d_view<real_type, DT>(do_not_init_tag("instance rates"), n_cells);
  _instance_rates_scan = value_type_1d_view<real_type, DT>(do_not_init_tag("instance rates"), n_cells + 1);

  /// set interaction range
  _n_cells_interaction_x = n_cells_interaction_x;
  _n_cells_interaction_y = n_cells_interaction_y;

  /// update rates
  Impl::updateEventRatesLatticeDevice(_sid, _lattice, _dictionary, _instance_rates, verbose);

  //std::ofstream sum_rates_ofs;
  //sum_rates_ofs.open("serial_test_rates.txt", std::ios::out | std::ios::trunc);
  //sum_rates_ofs << std::setprecision(15);
  //for (ordinal_type i = 0; i < _instance_rates.extent(0); i++) {
  //  sum_rates_ofs << real_type(_instance_rates(i)) << '\n';
  //}

  //sum_rates_ofs.close();

  //std::ofstream proc_count_ofs;
  //proc_count_ofs.open("serial_proc_counts.txt", std::ios::out | std::ios::trunc);
  //for (ordinal_type i = 0, iend = _dictionary._found_process_counter.extent(0); i < iend; i++) {
  //  proc_count_ofs << _dictionary._found_process_counter(i) << '\n';
  //}
  //proc_count_ofs.close();
  //std::ofstream inst_count_ofs;
  //inst_count_ofs.open("serial_inst_counts.txt", std::ios::out | std::ios::trunc);
  //for (ordinal_type i = 0, iend = _dictionary._found_instance_counter.extent(0); i < iend; i++) {
  //  inst_count_ofs << _dictionary._found_instance_counter(i) << '\n';
  //}
  //inst_count_ofs.close();

  return 0;
}

template <typename DT>
ordinal_type SolverSerialRTA<DT>::advance(const value_type_2d_view<real_type, DT> t_in, const real_type t_step,
                                          const ordinal_type n_kmc_steps, const ordinal_type verbose) {
  const FunctionScope func_scope("SolverSerialRTA::advance", __FILE__, __LINE__, verbose);

  KINCAT_CHECK_ERROR(_max_kmc_steps_per_kernel_launch < n_kmc_steps,
                     "Error: SolverSerialRTA::advance n kmc steps is larger than max kmc steps per kernel launch");

  const ordinal_type last = _instance_rates.extent(0);
  for (ordinal_type iter = 0; iter < n_kmc_steps; ++iter) {
    /// scan rates
    _instance_rates_scan(0) = 0;
    for (ordinal_type i = 1, iend = _instance_rates_scan.extent(0); i < iend; ++i)
      _instance_rates_scan(i) = _instance_rates_scan(i - 1) + _instance_rates(i - 1);

    if (verbose > 2) {
      std::cout << "scanned instance rates: \n";
      for (ordinal_type i = 0, iend = _instance_rates.extent(0); i < iend; ++i) {
        ordinal_type k0, k1;
        _lattice.getLatticeCellIndex(i, k0, k1);
        std::cout << "[ [" << k0 << "," << k1 << "], [" << _instance_rates_scan(i) << ", "
                  << _instance_rates_scan(i + 1) << "] ]\n";
      }
    }

    /// compute random number, sum of rates, and time increment
    const real_type pi = getRandomNumber(verbose);
    const real_type zeta = getRandomNumber(verbose);
    const real_type sum_rates = _instance_rates_scan(last);
    const real_type dt = -std::log(zeta) / sum_rates;
    if (verbose > 2) {
      std::cout << "t: " << t_in(0, 0) << ", pi: " << pi << ", zeta: " << zeta << ", sum rates: " << sum_rates
                << ", dt: " << dt << "\n";
    }

    if (sum_rates == real_type(0)) {
      if (verbose > 2) {
        std::cout << "break the loop; sum rates is zero\n";
      }
      break;
    }

    /// select cell location
    const real_type pos = sum_rates * pi;
    const auto it = std::lower_bound(&_instance_rates_scan[0], &_instance_rates_scan[last], pos);

    /// check cell rate range includes the random variable
    const ordinal_type cid = std::distance(&_instance_rates_scan[0], it) - 1;
    KINCAT_CHECK_ERROR(pos < _instance_rates_scan(cid) || pos > _instance_rates_scan(cid + 1),
                       "Error: the position value is not located in the range of scan(cid) and scan(cid+1)");

    /// lattice index
    ordinal_type k0, k1;
    _lattice.getLatticeCellIndex(cid, k0, k1);

    if (verbose > 2) {
      std::cout << "position: " << pos << ", cid: " << cid << ", lattice cell index: [" << k0 << ", " << k1 << "]\n";
    }

    /// find instance and pattern
    ordinal_type idx_instance(-1), idx_pattern(-1), idx_variant(-1);
    real_type sum_rates_cell(-1);
    Impl::findEvent(_sid, _lattice, _dictionary, k0, k1, pos - _instance_rates_scan[cid], idx_instance, idx_pattern,
                    idx_variant, sum_rates_cell, verbose);
    const auto idx_variant_output = _dictionary._constraints(idx_instance, idx_variant);

    if (verbose > 2) {
      std::cout << "instance index: " << idx_instance << ", pattern index: " << idx_pattern
                << ", input variant index:" << idx_variant << ", output variant index: " << idx_variant_output << "\n";
    }

    if (idx_variant_output < 0) {
      /// do nothing, no event occurs this step
    } else {
      /// check the data is consistent
      const real_type epsilon(1e-10);
      KINCAT_CHECK_ERROR(std::abs(sum_rates_cell - _instance_rates(cid)) > epsilon,
                         "Error: rate on this site does not match to instance rate data");
      Impl::updateLatticeConfiguration(_sid, _lattice, _dictionary, k0, k1, idx_instance, idx_pattern, idx_variant,
                                       verbose);
      /// update rates
      Impl::updateEventRates(_sid, _lattice, _dictionary, k0, k1, _n_cells_interaction_x, _n_cells_interaction_y,
                             _instance_rates, verbose);
      /// record counter
      _counter.update(_sid, 0, idx_instance); //Only one domain...

      /// advance global time
      t_in(0, 0) += dt;

      /// break if t is bigger than t step
      if (t_step > 0 && t_in(0, 0) > t_step)
        break;
    }
  }
  return 0;
}

template <typename DT>
std::ostream &SolverSerialRTA<DT>::showMe(std::ostream &os, const std::string &label,
                                          const ordinal_type verbose) const {
  os << "-- Solver : " << label << "\n";

  return os;
}

/// eti for individual device type
#if defined(KOKKOS_ENABLE_SERIAL)
template struct SolverSerialRTA<typename UseThisDevice<Kokkos::Serial>::type>;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template struct SolverSerialRTA<typename UseThisDevice<Kokkos::OpenMP>::type>;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template struct SolverSerialRTA<typename UseThisDevice<Kokkos::Cuda>::type>;
#endif

} // namespace KinCat