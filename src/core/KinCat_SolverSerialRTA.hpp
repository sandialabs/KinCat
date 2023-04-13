#ifndef __KINCAT_SOLVER_SERIAL_RTA_HPP__
#define __KINCAT_SOLVER_SERIAL_RTA_HPP__

// clang-format off
#include "KinCat_Util.hpp"
#include "KinCat_SolverBase.hpp"
// clang-format on

namespace KinCat {

template <typename DeviceType> struct SolverSerialRTA : public SolverBase<DeviceType> {
private:
  ordinal_type _sid;

  value_type_1d_view<real_type, DeviceType> _instance_rates;      /// size of lattice
  value_type_1d_view<real_type, DeviceType> _instance_rates_scan; /// binary search

  ordinal_type _iter_random_pool;
  Kokkos::Random_XorShift64_Pool<DeviceType> _random;
  value_type_1d_view<real_type, DeviceType> _random_pool;

  ordinal_type _n_cells_interaction_x, _n_cells_interaction_y;
  Lattice<DeviceType> _lattice;
  ProcessDictionary<DeviceType> _dictionary;
  ProcessCounter<DeviceType> _counter;

  ordinal_type _max_kmc_steps_per_kernel_launch;

  ordinal_type fillRandomPool(const ordinal_type verbose = 0);
  real_type getRandomNumber(const ordinal_type verbose = 0);

public:
  SolverSerialRTA();
  SolverSerialRTA(const SolverSerialRTA &b);
  virtual ~SolverSerialRTA();

  virtual ordinal_type initialize(const bool is_solver_random_number_variation, const ordinal_type random_seed,
                                  const ordinal_type random_pool_size,
                                  const ordinal_type max_kmc_steps_per_kernel_launch,
                                  const ordinal_type n_cells_interaction_x, const ordinal_type n_cells_interaction_y,
                                  const Lattice<DeviceType> &lattice, const ProcessDictionary<DeviceType> &dictionary,
                                  const ProcessCounter<DeviceType> &counter, const ordinal_type verbose = 0);
  virtual ordinal_type advance(const value_type_2d_view<real_type, DeviceType> t_in, const real_type t_step,
                               const ordinal_type n_kmc_steps, const ordinal_type verbose = 0);
  virtual ordinal_type advance(const value_type_1d_view<real_type, DeviceType> t_in, const real_type t_step,
                               const ordinal_type n_kmc_steps, const ordinal_type verbose = 0) {
    KINCAT_CHECK_ERROR(true, "Error: this method is for batch interface and not supported in SerialRTA");
  }
  virtual std::ostream &showMe(std::ostream &os, const std::string &label, const ordinal_type verbose = 0) const;
};

} // namespace KinCat

#endif
