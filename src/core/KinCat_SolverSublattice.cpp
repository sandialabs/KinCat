// clang-format off
#include "KinCat_Util.hpp"
#include "KinCat_SolverSublattice.hpp"
#include "KinCat_SolverInternal.hpp"
#include <chrono>
// clang-format on

namespace KinCat {

template <typename DT> SolverSublattice<DT>::SolverSublattice() = default;
template <typename DT> SolverSublattice<DT>::SolverSublattice(const SolverSublattice &b) = default;
template <typename DT> SolverSublattice<DT>::~SolverSublattice() = default;

template <typename DT> ordinal_type SolverSublattice<DT>::fillRandomPool(const ordinal_type verbose) {
  const FunctionScope func_scope("SolverSublattice::fillRandomPool", __FILE__, __LINE__, verbose);
  Kokkos::fill_random(_random_pool, _random, real_type(1));
  Kokkos::deep_copy(_iter_random_pool, 0);

  return 0;
}

template <typename DT>
ordinal_type SolverSublattice<DT>::initialize(const bool is_solver_random_number_variation,
                                              const ordinal_type random_seed, const ordinal_type random_pool_size,
                                              const ordinal_type max_kmc_steps_per_kernel_launch,
                                              const ordinal_type n_cells_interaction_x,
                                              const ordinal_type n_cells_interaction_y, const Lattice<DT> &arg_lattice,
                                              const ProcessDictionary<DT> &arg_dictionary,
                                              const ProcessCounter<DT> &arg_counter, const ordinal_type verbose) {
  const FunctionScope func_scope("SolverSublattice::initialize", __FILE__, __LINE__, verbose);

  _sid = 0;

  /// set member variables
  _lattice = arg_lattice;
  _dictionary = arg_dictionary;
  _counter = arg_counter;

  /// set random variable
  _random = Kokkos::Random_XorShift64_Pool<DT>(random_seed);

  /// set max kmc steps per kernel launch (used for workspace and random pool setup)
  _max_kmc_steps_per_kernel_launch = max_kmc_steps_per_kernel_launch;

  /// subdomain check
  ordinal_type n_domains, n_domains_x, n_domains_y;
  n_domains = _lattice.getNumberOfDomains(n_domains_x, n_domains_y);

  const ordinal_type n_cells_x = _lattice._n_cells_domain_x * n_domains_x;
  const ordinal_type n_cells_y = _lattice._n_cells_domain_y * n_domains_y;
  KINCAT_CHECK_ERROR(_lattice._n_cells_x != n_cells_x, "Error: latice length x is not multiple of subdomain x");
  KINCAT_CHECK_ERROR(_lattice._n_cells_y != n_cells_y, "Error: latice length y is not multiple of subdomain y");

  /// set interaction range
  _n_cells_interaction_x = n_cells_interaction_x;
  _n_cells_interaction_y = n_cells_interaction_y;
  KINCAT_CHECK_ERROR(_lattice._n_cells_x != n_cells_x, "Error: latice length x is not multiple of subdomain x");
  KINCAT_CHECK_ERROR(_lattice._n_cells_x != n_cells_x, "Error: latice length x is not multiple of subdomain x");

  /// uniform random values (we do not want to refill random numbers during a kmc run)
  const ordinal_type required_random_pool_size = _max_kmc_steps_per_kernel_launch * 2;
  // const ordinal_type p0 = random_pool_size / required_random_pool_size;
  // const ordinal_type p1 = random_pool_size % required_random_pool_size;

  /// ignore the random pool size argument (use the minimum per kernel launch)
  // const ordinal_type actual_random_pool_size = (p0 + (p1 > 0)) * required_random_pool_size;
  const ordinal_type actual_random_pool_size = required_random_pool_size;

  /// uniform random values
  _iter_random_pool =
      value_type_2d_view<ordinal_type, DT>(do_not_init_tag("iter random pool"), n_domains_x, n_domains_y);
  _random_pool = value_type_3d_view<real_type, DT>(do_not_init_tag("random pool"), n_domains_x, n_domains_y,
                                                   actual_random_pool_size);

  /// create instance list
  _instance_rates = value_type_3d_view<real_type, DT>(do_not_init_tag("instance rates"), n_domains_x, n_domains_y,
                                                      _lattice.getNumberOfDomainCells());
  _instance_rates_scan = value_type_3d_view<real_type, DT>(do_not_init_tag("instance rates"), n_domains_x, n_domains_y,
                                                           _instance_rates.extent(2) + 1);

  /// update rates
  Impl::updateEventRatesLatticeDevice(_sid, _lattice, _dictionary, _instance_rates, verbose);

  return 0;
}

template <typename DT>
ordinal_type SolverSublattice<DT>::advance(const value_type_2d_view<real_type, DT> t_in, const real_type t_step,
                                           const ordinal_type n_kmc_steps, const ordinal_type verbose) {
  const FunctionScope func_scope("SolverSublattice::advance", __FILE__, __LINE__, verbose);

  const ordinal_type quad[8] = {0, 0, 1, 0, 0, 1, 1, 1};

  using policy_type = Kokkos::TeamPolicy<typename DT::execution_space>;

  ordinal_type n_domains, n_domains_x, n_domains_y;
  n_domains = _lattice.getNumberOfDomains(n_domains_x, n_domains_y);
  const ordinal_type n_domains_x2 = n_domains_x / 2, n_domains_y2 = n_domains_y / 2,
                     n_domains_4 = n_domains_x2 * n_domains_y2;

  const policy_type policy(n_domains_4, Kokkos::AUTO(), Kokkos::AUTO());

  const ordinal_type sid = _sid;

  const auto lattice = _lattice;
  const auto dictionary = _dictionary;
  const auto counter = _counter;

  const auto instance_rates_all = _instance_rates;
  const auto instance_rates_scan_all = _instance_rates_scan;
  const auto iter_random_pool_all = _iter_random_pool;
  const auto random_pool_all = _random_pool;

  const ordinal_type instance_rates_extent = instance_rates_all.extent(2);
  const ordinal_type n_cells_interaction_x = _n_cells_interaction_x;
  const ordinal_type n_cells_interaction_y = _n_cells_interaction_y;

  KINCAT_CHECK_ERROR(fillRandomPool(verbose), "Error: fillRandomPool returns non-zero error code");

  unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 random_number_generator(seed1);

  /// loop over four quadrant sublattices, each with unique 'colors'.
  /// need to ensure order of quadrant looping changes.
  bool kernel_done = false;

  std::vector<ordinal_type> q_order = {0, 1, 2, 3}; // colors

  std::shuffle(q_order.begin(), q_order.end(), random_number_generator);
  // std::cout << "Quad Order: " << q_order[0] << q_order[1] << q_order[2] << q_order[3] << '\n';
  for (ordinal_type q = 0; q < 4; ++q) {              // loop over colors
    const ordinal_type q0 = quad[q_order[q] * 2 + 0]; // getting coordinates of quadrant
    const ordinal_type q1 = quad[q_order[q] * 2 + 1];
    /// launch a kernel for all of this color's subdomains
    for (ordinal_type iter = 0; iter < n_kmc_steps; ++iter) {
      ordinal_type n_domains_step(0);
      Kokkos::Sum<ordinal_type> reducer_value(n_domains_step);
      Kokkos::parallel_reduce(
          policy,
          KOKKOS_LAMBDA(const typename policy_type::member_type &member, ordinal_type &update) {
            const ordinal_type did = member.league_rank();
            const ordinal_type d0_2 = did / n_domains_y2, d1_2 = did % n_domains_y2;
            const ordinal_type d0 = 2 * d0_2 + q0, d1 = 2 * d1_2 + q1;
            Kokkos::single(Kokkos::PerTeam(member), [&]() { update += (t_step > 0 ? (t_in(d0, d1) >= t_step) : 0); });

            /// process the selected domain
            if (t_step == 0 || (t_step > 0 && t_in(d0, d1) < t_step)) {
              const value_type_1d_view<real_type, DT> instance_rates =
                  Kokkos::subview(instance_rates_all, d0, d1, Kokkos::ALL());
              const value_type_1d_view<real_type, DT> instance_rates_scan =
                  Kokkos::subview(instance_rates_scan_all, d0, d1, Kokkos::ALL());

              ordinal_type &iter = iter_random_pool_all(d0, d1);
              const auto random_pool = Kokkos::subview(random_pool_all, d0, d1, Kokkos::ALL());

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
                const ordinal_type lid = static_cast<ordinal_type>(loc - first) - 1;

                ordinal_type l0, l1;
                lattice.getDomainLatticeCellIndex(lid, l0, l1);

                ordinal_type k0, k1;
                lattice.getLatticeCellIndex(d0, d1, l0, l1, k0, k1);

                /// find instance and pattern
                ordinal_type idx_instance(-1), idx_pattern(-1), idx_variant(-1);
                real_type sum_rates_cell(-1);
                Impl::findEvent(sid, lattice, dictionary, k0, k1, pos - instance_rates_scan(lid), idx_instance,
                                idx_pattern, idx_variant, sum_rates_cell, verbose);

                /// check instance is possible or not
                const auto idx_variant_output = dictionary._constraints(idx_instance, idx_variant);

                if (idx_variant_output < 0) {
                  /// do nothing, no event occurs this step
                } else if (((t_in(d0, d1) + dt) > t_step) && (t_step > 0.0)) {
                  /// event would occur beyond t_step and thus is rejected, but local time is still evolved up to t_step
                  t_in(d0, d1) = t_step;
                } else {
                  /// check the data is consistent
                  const real_type epsilon(1e-10);
                  KINCAT_CHECK_ERROR(ats<real_type>::abs(sum_rates_cell - instance_rates(lid)) > epsilon,
                                     "Error: rate on this site does not match to instance rate data");

                  /// update lattice configuration
                  Impl::updateLatticeConfiguration(sid, lattice, dictionary, k0, k1, idx_instance, idx_pattern,
                                                   idx_variant, verbose);
                  Impl::updateEventRates(sid, lattice, dictionary, k0, k1, n_cells_interaction_x, n_cells_interaction_y,
                                         instance_rates_all, verbose);

                  /// update counter
                  counter.update(sid, idx_instance);

                  /// advance time
                  t_in(d0, d1) += dt;
                }
              });
            }
          },
          reducer_value);
      typename DT::execution_space().fence(); // TODO: Check if necessary
      /// check if all subdomains for this color are updated up to t_step -- if so, we are done with the color
      if (n_domains_step == n_domains_4) {
        break;
      }
    }
  }

  return 0;
}

template <typename DT>
std::ostream &SolverSublattice<DT>::showMe(std::ostream &os, const std::string &label,
                                           const ordinal_type verbose) const {
  os << "-- Solver : " << label << "\n";

  return os;
}

/// eti for individual device type
#if defined(KOKKOS_ENABLE_SERIAL)
template struct SolverSublattice<typename UseThisDevice<Kokkos::Serial>::type>;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template struct SolverSublattice<typename UseThisDevice<Kokkos::OpenMP>::type>;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template struct SolverSublattice<typename UseThisDevice<Kokkos::Cuda>::type>;
#endif

} // namespace KinCat
