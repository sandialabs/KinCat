// clang-format off
#include "KinCat_Util.hpp"
#include "KinCat_SolverBatchRTA.hpp"
#include "KinCat_SolverInternal.hpp"
// clang-format on

namespace KinCat {

template <typename DT> SolverBatchRTA<DT>::SolverBatchRTA() = default;
template <typename DT> SolverBatchRTA<DT>::SolverBatchRTA(const SolverBatchRTA &b) = default;
template <typename DT> SolverBatchRTA<DT>::~SolverBatchRTA() = default;

template <typename DT> ordinal_type SolverBatchRTA<DT>::fillRandomPool(const ordinal_type verbose) {
  const FunctionScope func_scope("SolverBatchRTA::fillRandomPool", __FILE__, __LINE__, verbose);
  Kokkos::fill_random(_random_pool, _random, real_type(1));
  Kokkos::deep_copy(_iter_random_pool, 0);

  return 0;
}

template <typename DT>
ordinal_type
SolverBatchRTA<DT>::initialize(const bool is_solver_random_number_variation, const ordinal_type random_seed,
                               const ordinal_type random_pool_size, const ordinal_type max_kmc_steps_per_kernel_launch,
                               const ordinal_type n_cells_interaction_x, const ordinal_type n_cells_interaction_y,
                               const Lattice<DT> &arg_lattice, const ProcessDictionary<DT> &arg_dictionary,
                               const ProcessCounter<DT> &arg_counter, const ordinal_type verbose) {
  const FunctionScope func_scope("SolverBatchRTA::initialize", __FILE__, __LINE__, verbose);

  const ordinal_type n_samples = arg_lattice.getNumberOfSamples();

  /// set member variables
  _lattice = arg_lattice;
  _dictionary = arg_dictionary;
  _counter = arg_counter;

  /// serial version always has the same domain as the lattice
  _lattice.setDomain(_lattice._n_cells_x, _lattice._n_cells_y);

  /// set random variable
  _random = Kokkos::Random_XorShift64_Pool<DT>(random_seed);

  /// set max kmc steps per kernel launch (used for workspace and random pool setup)
  _max_kmc_steps_per_kernel_launch = max_kmc_steps_per_kernel_launch;

  /// uniform random values (we do not want to refill random numbers during a kmc run)
  const ordinal_type required_random_pool_size = _max_kmc_steps_per_kernel_launch * 2;
  const ordinal_type actual_random_pool_size = required_random_pool_size;

  _iter_random_pool = value_type_1d_view<ordinal_type, DT>(do_not_init_tag("iter random pool"), n_samples);

  _random_pool = value_type_2d_view<real_type, DT>(
      do_not_init_tag("random pool"), (is_solver_random_number_variation ? n_samples : 1), actual_random_pool_size);
  KINCAT_CHECK_ERROR(fillRandomPool(verbose), "Error: fillRandomPool returns non-zero error code");

  /// create process instance list
  const ordinal_type n_cells = _lattice.getNumberOfCells();
  _instance_rates = value_type_2d_view<real_type, DT>(do_not_init_tag("instance rates"), n_samples, n_cells);
  _instance_rates_scan = value_type_2d_view<real_type, DT>(do_not_init_tag("instance rates"), n_samples, n_cells + 1);

  /// set interaction range
  _n_cells_interaction_x = n_cells_interaction_x;
  _n_cells_interaction_y = n_cells_interaction_y;

  /// update rates
  Impl::updateEventRatesLatticeDevice(_lattice, _dictionary, _instance_rates, verbose);

  return 0;
}

template <typename DT>
ordinal_type SolverBatchRTA<DT>::advance(const value_type_1d_view<real_type, DT> t_in, const real_type t_step,
                                         const ordinal_type n_kmc_steps, const ordinal_type verbose) {
  const FunctionScope func_scope("SolverBatchRTA::advance", __FILE__, __LINE__, verbose);

  using policy_type = Kokkos::TeamPolicy<typename DT::execution_space>;

  const ordinal_type n_samples = _lattice.getNumberOfSamples();
  const policy_type policy(n_samples, Kokkos::AUTO(), Kokkos::AUTO());

  const auto lattice = _lattice;
  const auto dictionary = _dictionary;
  const auto counter = _counter;

  const auto instance_rates_all = _instance_rates;
  const auto instance_rates_scan_all = _instance_rates_scan;
  const auto iter_random_pool_all = _iter_random_pool;
  const auto random_pool_all = _random_pool;

  const ordinal_type instance_rates_extent = instance_rates_all.extent(1);
  const ordinal_type n_cells_interaction_x = _n_cells_interaction_x;
  const ordinal_type n_cells_interaction_y = _n_cells_interaction_y;

  KINCAT_CHECK_ERROR(fillRandomPool(verbose), "Error: fillRandomPool returns non-zero error code");

  /// TODO:: consider to put outer loop inside of parallel region so that we can launch
  ///        a single kernel running multiple kmc steps
  for (ordinal_type iter = 0; iter < n_kmc_steps; ++iter) {
    /// launch a kernel for all samples
    ordinal_type n_samples_step(0);
    Kokkos::Sum<ordinal_type> reducer_value(n_samples_step);
    Kokkos::parallel_reduce(
        policy,
        KOKKOS_LAMBDA(const typename policy_type::member_type &member, ordinal_type &update) {
          const ordinal_type sid = member.league_rank();
          Kokkos::single(Kokkos::PerTeam(member), [&]() { update += (t_step > 0 ? (t_in(sid) >= t_step) : 0); });

          /// if t_step is zero, then the condition to reach t_step is ignored
          /// simulation runs upto the n kmc steps
          if (t_step == 0 || (t_step > 0 && t_in(sid) < t_step)) {
            /// process the selected sample
            const value_type_1d_view<real_type, DT> instance_rates =
                Kokkos::subview(instance_rates_all, sid, Kokkos::ALL());
            const value_type_1d_view<real_type, DT> instance_rates_scan =
                Kokkos::subview(instance_rates_scan_all, sid, Kokkos::ALL());

            ordinal_type &iter = iter_random_pool_all(sid);
            const ordinal_type random_pool_n_samples = random_pool_all.extent(0);
            const auto random_pool = Kokkos::subview(random_pool_all, sid % random_pool_n_samples, Kokkos::ALL());

            /// exclusive scan
            Kokkos::parallel_scan(Kokkos::ThreadVectorRange(member, instance_rates_extent + 1),
                                  [=](const ordinal_type &l, real_type &update, bool final) {
                                    if (final)
                                      instance_rates_scan(l) = update;

                                    if (l < instance_rates_extent)
                                      update += instance_rates(l);
                                  });
            member.team_barrier();

            Kokkos::single(Kokkos::PerTeam(member), [&]() {
              /// compute random number, sum of rates, and time increment
              const real_type pi = random_pool(iter++);
              const real_type zeta = random_pool(iter++);
              const real_type sum_rates = instance_rates_scan(instance_rates_extent);
              const real_type dt = -ats<real_type>::log(zeta) / sum_rates;

              KINCAT_CHECK_ERROR(sum_rates == real_type(0), "Error: sum rates is zero");

              const real_type pos = sum_rates * pi;
              real_type *first = &instance_rates_scan[0];
              real_type *last = &instance_rates_scan[instance_rates_extent];
              real_type *loc =
                  Impl::lower_bound(first, last, pos, [](real_type left, real_type right) { return left < right; });
              /// check cell rate range includes the random variable
              const ordinal_type cid = static_cast<ordinal_type>(loc - first) - 1;

              ordinal_type k0, k1;
              lattice.getLatticeCellIndex(cid, k0, k1);

              /// find instance and pattern
              ordinal_type idx_instance(-1), idx_pattern(-1), idx_variant(-1);
              real_type sum_rates_cell(-1);
              Impl::findEvent(sid, lattice, dictionary, k0, k1, pos - instance_rates_scan(cid), idx_instance,
                              idx_pattern, idx_variant, sum_rates_cell, verbose);

              /// check instance is possible or not
              const auto idx_variant_output = dictionary._constraints(idx_instance, idx_variant);

              if (idx_variant_output < 0) {
                /// do nothing this step
              } else {
                /// check the data is consistent
                const real_type epsilon(1e-10);
                KINCAT_CHECK_ERROR(ats<real_type>::abs(sum_rates_cell - instance_rates(cid)) > epsilon,
                                   "Error: rate on this site does not match to instance rate data");

                /// update lattice configuration
                Impl::updateLatticeConfiguration(sid, lattice, dictionary, k0, k1, idx_instance, idx_pattern,
                                                 idx_variant, verbose);
                Impl::updateEventRates(sid, lattice, dictionary, k0, k1, n_cells_interaction_x, n_cells_interaction_y,
                                       instance_rates, verbose);

                /// update counter
                counter.update(sid, 0, idx_instance); //currently only uses serial algorithm, so only one domain

                /// advance time
                t_in(sid) += dt;
              }
            });
          }
        },
        reducer_value);

    /// all samples reach to t_step
    if (n_samples_step == n_samples)
      break;
  }

  return 0;
}

template <typename DT>
std::ostream &SolverBatchRTA<DT>::showMe(std::ostream &os, const std::string &label, const ordinal_type verbose) const {
  os << "-- Solver : " << label << "\n";
  return os;
}

/// eti for individual device type
#if defined(KOKKOS_ENABLE_SERIAL)
template struct SolverBatchRTA<typename UseThisDevice<Kokkos::Serial>::type>;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template struct SolverBatchRTA<typename UseThisDevice<Kokkos::OpenMP>::type>;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template struct SolverBatchRTA<typename UseThisDevice<Kokkos::Cuda>::type>;
#endif

} // namespace KinCat
