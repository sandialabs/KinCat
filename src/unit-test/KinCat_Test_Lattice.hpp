#ifndef __KINCAT_TEST_LATTICE_HPP__
#define __KINCAT_TEST_LATTICE_HPP__

#include "KinCat_Lattice.hpp"
#include "KinCat_Main_Util.hpp"

TEST(Lattice, CommonInput) {
  {
    const std::vector<std::string> input = {"test-files/input.json", "test-files/input-batch.json"};
    for (auto &filename : input) {
      std::cerr << ">> Testing filename : " << filename << "\n";
      KinCat::Input in;
      in.parse(filename);

      /// Test various lattices?///
      // Input
      KinCat::Lattice<device_type> lattice(1, in._n_cells_x, in._n_cells_y, in._edge_vectors, in._basis_vectors,
                                           in._n_species);
      EXPECT_EQ(lattice.getNumberOfSites(), lattice._sites.extent(1));
    }
  }
}

#endif
