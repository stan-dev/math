#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <stdexcept>
#include <vector>

// Just test the use case where a single matrix is passed in to VectorViewMvt.
// The case where a std::vector of matrices is passed in is broken and will be
// refactored.

TEST(MetaTraits, VectorViewMvt_vector_double) {
  using stan::VectorViewMvt;
  using Eigen::VectorXd;
  using Eigen::RowVectorXd;
  using Eigen::Dynamic;

  VectorXd a(10);
  for (size_t n = 0; n < 10; ++n)
    a[n] = n;
  VectorViewMvt<VectorXd> av(a);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(a[n], av[0][n]);
  for (size_t n = 0; n < 10; ++n)
    av[120*n][n] = n + 10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10 + n, av[0][n]);
    EXPECT_FLOAT_EQ(10 + n, a[n]);
  }

  const VectorXd b(a);
  VectorViewMvt<const VectorXd> bv(b);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(b[n], bv[0][n]);

  RowVectorXd c(10);
  for (size_t n = 0; n < 10; ++n)
    c[n] = n;
  VectorViewMvt<RowVectorXd> cv(c);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(c[n], cv[0][n]);
  for (size_t n = 0; n < 10; ++n)
    cv[0][n] = n + 10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10 + n, cv[0][n]);
    EXPECT_FLOAT_EQ(10 + n, c[n]);
  }

  const RowVectorXd d(c);
  VectorViewMvt<const RowVectorXd, stan::is_vector<RowVectorXd>::value > dv(d);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(d[n], dv[0][n]);
}

TEST(MetaTraits, VectorViewMvt_matrix_double) {
  using stan::VectorViewMvt;
  using Eigen::MatrixXd;
  using Eigen::Dynamic;

  MatrixXd a(10, 10);
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 10; ++j) {
      a(i, j) = 100 * i + j;
    }
  }
  VectorViewMvt<MatrixXd> av(a);
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 10; ++j) {
      EXPECT_FLOAT_EQ(av[0](i, j), a(i, j));
    }
  }

  //Test mutation passes thru to original matrix
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 10; ++j) {
      av[0](i, j) = a(i, j) + 10;
    }
  }
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 10; ++j) {
      EXPECT_FLOAT_EQ(av[0](i, j), a(i, j));
    }
  }

  const MatrixXd b(a);
  VectorViewMvt<const MatrixXd> bv(b);
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 10; ++j) {
      EXPECT_FLOAT_EQ(bv[0](i, j), b(i, j));
    }
  }
}
