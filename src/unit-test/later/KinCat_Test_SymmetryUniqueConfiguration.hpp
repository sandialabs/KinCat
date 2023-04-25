#ifndef __KINCAT_TEST_SYMMETRY_UNIQUE_CONFIGURATION_HPP__
#define __KINCAT_TEST_SYMMETRY_UNIQUE_CONFIGURATION_HPP__

TEST(SymmetryUniqueConfiguration, ConstructorAndDestructor) {
  {
    auto sym_unique_conf = new KinCat::SymmetryUniqueConfiguration();
    delete sym_unique_conf;
  }
}

TEST(SymmetryUniqueConfiguration, buildUnitaryMatrices) {
  {
    using Trans = KinCat::Trans;
    using ats = KinCat::ArithTraits<real_type>;

    KinCat::SymmetryUniqueConfiguration sym_unique_conf;

    const ordinal_type dim(2);
    sym_unique_conf.buildUnitaryMatrices(dim);
    const auto unitary_matrices = sym_unique_conf.getUnitaryMatrices();

    const ordinal_type n = unitary_matrices.extent(0);
    EXPECT_EQ(n, 24);
    // printf(" # of unique unitary matrices %d\n", n);

    const real_type eps = 10 * ats::epsilon();
    for (ordinal_type k = 0; k < n; ++k) {
      const auto &A = unitary_matrices(k);

      /// unitary matrix
      const real_type det = ats::abs(determinant(A)), one(1), zero(0);
      EXPECT_NEAR(det, one, eps);

      /// unitary matrix recovers the identity
      bool is_identity(false);
      KinCat::real_type_matrix_3d Ac(A);
      for (ordinal_type i = 0; i < 6; ++i) {
        multiply(one, Trans::NoTranspose(), Ac, Trans::NoTranspose(), A, zero, Ac);
        if (KinCat::isIdentity(Ac)) {
          is_identity = true;
          break;
        }
      }
      EXPECT_TRUE(is_identity);
    }
  }
}

TEST(SymmetryUniqueConfiguration, findSymmetryOperations) {
  {
    using ats = KinCat::ArithTraits<real_type>;
    using KinCat::Trans;

    KinCat::SymmetryUniqueConfiguration sym_unique_conf;

    /// testing prime system
    KinCat::real_type_matrix_3d C, C_inv;
    C(0, 0) = 3.0 / 2.0;
    C(1, 0) = std::sqrt(3) / 2.0;

    C(0, 1) = 0, 0;
    C(1, 1) = std::sqrt(3);

    C(2, 2) = 1.0;

    KinCat::inverse(C, C_inv);

    /// empty basis list
    Kokkos::View<KinCat::real_type_vector_3d *, KinCat::host_device_type> basis_list;

    sym_unique_conf.findSymmetryOperations(C, basis_list);
    const auto sym_op_array = sym_unique_conf.getSymmetryOperations();

    const ordinal_type n = sym_op_array.extent(0);
    EXPECT_EQ(n, 12);
    // printf(" # of symmetry operations %d\n", n);

    const real_type eps = 10 * ats::epsilon();
    for (ordinal_type k = 0; k < n; ++k) {
      const auto &sym = sym_op_array(k);
      const auto &M = sym.getRotation();

      KinCat::real_type_matrix_3d A;

      /// A := C*M*C_inv;
      const real_type one(1), zero(0);

      KinCat::multiply(one, Trans::NoTranspose(), M, Trans::NoTranspose(), C_inv, zero, A);
      KinCat::multiply(one, Trans::NoTranspose(), C, Trans::NoTranspose(), A, zero, A);

      /// A := A*A^T
      KinCat::multiply(one, Trans::NoTranspose(), A, Trans::Transpose(), A, zero, A);
      const bool is_identity = KinCat::isIdentity(A);

      EXPECT_TRUE(is_identity);

      /// TODO: test basis list
    }
  }
}

#endif
