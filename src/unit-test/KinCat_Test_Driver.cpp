#include "KinCat.hpp"
#include <gtest/gtest.h>

using real_type = KinCat::real_type;
using ordinal_type = KinCat::ordinal_type;
using size_type = KinCat::size_type;
using site_type = KinCat::site_type;
using device_type = KinCat::device_type;
using host_device_type = KinCat::host_device_type;

#include "KinCat_Test_JSON.hpp"
#include "KinCat_Test_Lattice.hpp"
#include "KinCat_Test_Runs.hpp"

int main(int argc, char *argv[]) {
  int r_val(0);
  Kokkos::initialize(argc, argv);
  {
    const bool detail = false;
    KinCat::exec_space().print_configuration(std::cout, detail);
    KinCat::host_exec_space().print_configuration(std::cout, detail);

    ::testing::InitGoogleTest(&argc, argv);
    r_val = RUN_ALL_TESTS();
  }
  Kokkos::finalize();

  return r_val;
}
