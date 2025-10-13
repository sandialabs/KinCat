#ifndef __KINCAT_UTIL_HPP__
#define __KINCAT_UTIL_HPP__

/// basic std headers
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <cassert>
#include <cmath>
#include <ctime>
#include <limits>

#include <complex>
//#include <stdio.h>

/// kokkos
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "Kokkos_Timer.hpp"
//#include "Kokkos_UnorderedMap.hpp"

/// cmake auto-generated config file
#include "KinCat_Config.hpp"

namespace KinCat {
/// global variable types
using real_type = double;
using ordinal_type = int;
using size_type = size_t;

/// site data type
#if defined(KINCAT_ENABLE_SITE_TYPE_CHAR)
//#warning "site_type defined char/short"
using site_type = short; 
/// char provides up to 256 species using 1 byte, far more than needed. 
/// However, using char means that it may be compiled as either signed or unsigned char. 
/// Specifying either signed or unsigned leads to a compile error in Kokkos random.
/// May be able to use char if switch kincatpy to use large number e.g. 100, as empty site, rather than -1.
/// Choose instead to use short, which is a signed type. 
/// Currently using int for HDF5 compatibility. Need to write function for shorts. 
#elif defined(KINCAT_ENABLE_SITE_TYPE_SHORT)
//#warning "site_type defined short"
using site_type = short; /// using 2 byte //KK set to short int, CD changed
#else
//#warning "site_type defined int"
using site_type = int; /// using 4 byte, this is probably overkill
#endif

/// kokkos execution and device type
using exec_space = Kokkos::DefaultExecutionSpace;
using host_exec_space = Kokkos::DefaultHostExecutionSpace;

/// layout
using layout_right = Kokkos::LayoutRight;
using layout_left = Kokkos::LayoutLeft;

/// default device type
template <typename SpT> struct UseThisDevice;

#if defined(KOKKOS_ENABLE_CUDA)
template <> struct UseThisDevice<Kokkos::Cuda> {
  using type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
  using device_type = type;
};
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template <> struct UseThisDevice<Kokkos::OpenMP> {
  using type = Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>;
  using device_type = type;
};
#endif
#if defined(KOKKOS_ENABLE_SERIAL)
template <> struct UseThisDevice<Kokkos::Serial> {
  using type = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
  using device_type = type;
};
#endif

/// device types
using device_type = typename UseThisDevice<exec_space>::type;
using host_device_type = typename UseThisDevice<host_exec_space>::type;

using do_not_init_tag = Kokkos::ViewAllocateWithoutInitializing;

/// view types
template <typename ValueType, typename DeviceType>
using value_type_0d_view = Kokkos::View<ValueType, Kokkos::LayoutRight, DeviceType>;

template <typename ValueType, typename DeviceType>
using value_type_1d_view = Kokkos::View<ValueType *, Kokkos::LayoutRight, DeviceType>;

template <typename ValueType, typename DeviceType>
using value_type_2d_view = Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType>;

template <typename ValueType, typename DeviceType>
using value_type_3d_view = Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType>;

template <typename ValueType, typename DeviceType>
using value_type_4d_view = Kokkos::View<ValueType ****, Kokkos::LayoutRight, DeviceType>;

template <typename ValueType, typename DeviceType>
using value_type_5d_view = Kokkos::View<ValueType *****, Kokkos::LayoutRight, DeviceType>;

#if defined(KINCAT_ENABLE_DEBUG)
/// do actual error check; if the err is true, then it throw an exception
#if defined(__CUDA_ARCH__)
#define KINCAT_DEBUG_CHECK_ERROR(err, msg)                                                                             \
  if (err)                                                                                                             \
    Kokkos::abort(msg);
#else
#define KINCAT_DEBUG_CHECK_ERROR(err, msg)                                                                             \
  if (err)                                                                                                             \
    throw std::logic_error(msg);
#endif
#else
#define KINCAT_DEBUG_CHECK_ERROR(err, msg)
#endif

#if defined(__CUDA_ARCH__)
#define KINCAT_CHECK_ERROR(err, msg)                                                                                   \
  if (err)                                                                                                             \
    Kokkos::abort(msg);
#else
#define KINCAT_CHECK_ERROR(err, msg)                                                                                   \
  if (err)                                                                                                             \
    throw std::logic_error(msg);
#endif

template <typename T> struct ArithTraits;

template <> struct ArithTraits<float> {
  using value_type = float;
  using magnitude_type = float;
  static KOKKOS_FORCEINLINE_FUNCTION value_type log(const value_type &x) { return ::log(x); }
  static KOKKOS_FORCEINLINE_FUNCTION magnitude_type abs(const value_type &x) { return ::fabs(x); }
  static KOKKOS_FORCEINLINE_FUNCTION magnitude_type epsilon() { return DBL_EPSILON; }
};

template <> struct ArithTraits<double> {
  using value_type = double;
  using magnitude_type = double;
  static KOKKOS_FORCEINLINE_FUNCTION value_type log(const value_type &x) { return ::log(x); }
  static KOKKOS_FORCEINLINE_FUNCTION magnitude_type abs(const value_type &x) { return ::fabs(x); }
  static KOKKOS_FORCEINLINE_FUNCTION magnitude_type epsilon() { return DBL_EPSILON; }
};

template <> struct ArithTraits<int> {
  using value_type = int;
  using magnitude_type = int;
  using rats = ArithTraits<double>;
  static KOKKOS_FORCEINLINE_FUNCTION value_type log(const value_type &x) {
    return static_cast<value_type>(rats::log(x));
  }
  static KOKKOS_FORCEINLINE_FUNCTION magnitude_type abs(const value_type &x) { return x >= 0 ? x : -x; }
  static KOKKOS_FORCEINLINE_FUNCTION magnitude_type epsilon() { return 0; }
};

template <typename T> using ats = ArithTraits<T>;

struct FunctionScope {
  std::string _name;
  ordinal_type _verbose;
  FunctionScope(const std::string name, const std::string file, const ordinal_type lineno,
                const ordinal_type verbose = 0)
      : _name(name), _verbose(verbose) {
    if (_verbose > 1)
      std::cout << _name << " : begin";
    if (_verbose > 2)
      std::cout << ", file : " << file << ", line : " << lineno;
    if (_verbose > 1)
      std::cout << "\n";
  }
  ~FunctionScope() {
    if (_verbose > 1)
      std::cout << _name << " : end\n";
  }
#if defined(KINCAT_ENABLE_DEBUG)
  /// create actual stack operation for debugging purpose
#endif
};
} // namespace KinCat

#endif