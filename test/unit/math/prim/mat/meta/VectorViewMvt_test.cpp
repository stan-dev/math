#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <stdexcept>
#include <vector>

// Just test the use case where a single matrix is passed in to VectorViewMvt.
// The case where a std::vector of matrices is passed in is broken and will be refactored.

TEST(MetaTraits, VectorViewMvt_matrix_double) {
  using stan::VectorViewMvt;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<double,Dynamic,1> a(10);
  for (size_t n = 0; n < 10; ++n)
    a[n] = n;
  VectorViewMvt<Matrix<double,Dynamic,1> > av(a);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(a[n], av[0][n]);
  for (size_t n = 0; n < 10; ++n)
    av[120*n][n] = n+10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10+n, av[0][n]);
    EXPECT_FLOAT_EQ(10+n, a[n]);
  }

  const Matrix<double,Dynamic,1> b(a);
  VectorViewMvt<const Matrix<double,Dynamic,1> > bv(b);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(b[n], bv[0][n]);

  Matrix<double,1,Dynamic> c(10);
  for (size_t n = 0; n < 10; ++n)
    c[n] = n;
  VectorViewMvt<Matrix<double,1,Dynamic> > cv(c);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(c[n], cv[0][n]);
  for (size_t n = 0; n < 10; ++n)
    cv[0][n] = n+10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10+n, cv[0][n]);
    EXPECT_FLOAT_EQ(10+n, c[n]);
  }

  const Matrix<double,1,Dynamic> d(c);
  VectorViewMvt<const Matrix<double,1,Dynamic>,
          stan::is_vector<Matrix<double,1,Dynamic> >::value > dv(d);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(d[n], dv[0][n]);
}
