// clang-format off
#include "KinCat_SolverBase.hpp"
#include "KinCat_SolverSerialRTA.hpp"
#include "KinCat_SolverSublattice.hpp"
#include "KinCat_SolverBatchRTA.hpp"

// clang-format on

namespace KinCat {

template <typename DT>
SolverBase<DT> *SolverFactory<DT>::createSolver(const std::string &solver_type, const ordinal_type verbose) {
  const FunctionScope func_scope("SolverFactory<DT>::createSolver", __FILE__, __LINE__, verbose);

  SolverBase<DT> *r_val = nullptr;
  if (solver_type == "serial-rta") {
    r_val = new SolverSerialRTA<DT>();
  } else if (solver_type == "sublattice") {
    r_val = new SolverSublattice<DT>();
  } else if (solver_type == "batch-rta") {
    r_val = new SolverBatchRTA<DT>();
  } else if (solver_type == "timewarp") {
    // r_val = new SolverTimeWarp<DT>();
    KINCAT_CHECK_ERROR(true, "Error: TimeWarp is not implemented yet");
  }

  {
    std::stringstream ss;
    ss << "Error: failed to create a solver [" << solver_type << "]";
    KINCAT_CHECK_ERROR(r_val == nullptr, ss.str().c_str());
  }

  return r_val;
}

/// eti for individual device type
#if defined(KOKKOS_ENABLE_SERIAL)
template struct SolverFactory<typename UseThisDevice<Kokkos::Serial>::type>;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template struct SolverFactory<typename UseThisDevice<Kokkos::OpenMP>::type>;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template struct SolverFactory<typename UseThisDevice<Kokkos::Cuda>::type>;
#endif

} // namespace KinCat