#ifndef __KINCAT_TEST_JSON_HPP__
#define __KINCAT_TEST_JSON_HPP__

#include "KinCat_Main_Util.hpp"

TEST(JSON, CommonInput) {
  {
    const std::vector<std::string> input = {"test-files/input.json"}; //, "test-files/input-batch.json"};
    for (auto &filename : input) {
      std::cerr << ">> Testing filename : " << filename << "\n";
      KinCat::Input in;
      in.parse(filename, 4);

      /// event array match to rates array size
      EXPECT_EQ(in._processints.extent(0), in._rates.extent(1));

      /// number of samples from input rates should be 1
      EXPECT_EQ(in._rates.extent(0), 1);

      /// process index size should match to processes
      EXPECT_EQ(in._processes.size(), in._process_index.size());
      {
        /// check if process index map is correct
        for (ordinal_type i = 0, iend = in._processes.size(); i < iend; ++i) {
          EXPECT_EQ(in._process_index[in._processes[i]], i);
        }
      }
    }
  }
}

TEST(JSON, BatchSpecificInput) {
  {
    const std::string input_filename("test-files/input-batch.json");

    KinCat::Input in;
    in.parse(input_filename, 4);

    /// n samples should be set in batch input file
    std::cout << in._n_samples << "found samples \n";
    EXPECT_GT(in._n_samples, 0);

    /// solver should be batch-rta for now
    std::cout << in._solver_type << '\n';
    EXPECT_EQ(in._solver_type, "batch-rta");

    /// check ensemble input variation options
  }
}

TEST(JSON, EnsembleSpecificInput) {
  {
    // FINISH!!! TODO
    // Same as Batch Specific Input?
  }
}

#endif
