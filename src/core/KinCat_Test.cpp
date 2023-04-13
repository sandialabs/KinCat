#include "KinCat.hpp"
#include "KinCat_Util.hpp"

int main(int argc, char *argv[]) {
  using real_type = KinCat::real_type;
  using ordinal_type = KinCat::ordinal_type;

  Kokkos::initialize(argc, argv);
  {
    /// Craig's test
    printf("herp\n");
    const ordinal_type n_sites = 100;
    const real_type fill = 0.5;
    const ordinal_type n_fill_sites = fill * n_sites;

    // value_type_1d_view<ordinal_type> site_locations("site locations", n_fill_sites);
    // Kokkos::fill_random(site_locations, random, n_sites);
  }

  Kokkos::finalize();

  return 0;
}
