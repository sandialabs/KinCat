#ifndef __KINCAT_SOLVER_BASE_HPP__
#define __KINCAT_SOLVER_BASE_HPP__

// clang-format off
#include "KinCat_Util.hpp"
#include "KinCat_Lattice.hpp"
#include "KinCat_ProcessCounter.hpp"
#include "KinCat_ProcessDictionary.hpp"
// clang-format on

namespace KinCat {

template <typename DeviceType> struct SolverBase {
public:
  SolverBase() = default;
  SolverBase(const SolverBase &b) = default;
  virtual ~SolverBase() = default;

  virtual ordinal_type initialize(const bool is_solver_random_number_variation, const ordinal_type solver_random_seed,
                                  const ordinal_type solver_random_pool_size,
                                  const ordinal_type max_kmc_iterations_per_kernel_launch,
                                  const ordinal_type n_cells_interaction_x, const ordinal_type n_cells_interaction_y,
                                  const Lattice<DeviceType> &lattice, const ProcessDictionary<DeviceType> &dictionary,
                                  const ProcessCounter<DeviceType> &counter, const ordinal_type verbose = 0) = 0;

  /// single problem e.g., serial-rta, sublattice, timewarp
  virtual ordinal_type advance(const value_type_2d_view<real_type, DeviceType> t_in, const real_type t_step,
                               const ordinal_type n_kmc_iterations, const ordinal_type verbose = 0) = 0;

  /// batch problem e.g., batch-rta
  virtual ordinal_type advance(const value_type_1d_view<real_type, DeviceType> t_in, const real_type t_step,
                               const ordinal_type n_kmc_iterations, const ordinal_type verbose = 0) = 0;

  /// verbose show
  virtual std::ostream &showMe(std::ostream &os, const std::string &label, const ordinal_type verbose = 0) const = 0;
};

/// create a solver object and return the pointer
template <typename DeviceType> struct SolverFactory {
  static SolverBase<DeviceType> *createSolver(const std::string &solver_type, const ordinal_type verbose = 0);
};

} // namespace KinCat

#endif