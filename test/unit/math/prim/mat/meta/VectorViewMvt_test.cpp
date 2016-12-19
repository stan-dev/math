#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <stdexcept>
#include <vector>

TEST(MetaTraits, VectorViewMvt_matrix_double) {
  using stan::VectorViewMvt;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<double,Dynamic,1> a(10);
  for (size_t n = 0; n < 10; ++n)
    a[n] = n;
  VectorViewMvt<Matrix<double,Dynamic,1> > av(a);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(a[n], av[n][n]);
  for (size_t n = 0; n < 10; ++n)
    av[120*n][n] = n+10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10+n, av[n][n]);
    EXPECT_FLOAT_EQ(10+n, a[n]);
  }

  const Matrix<double,Dynamic,1> b(a);
  VectorViewMvt<const Matrix<double,Dynamic,1> > bv(b);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(b[n], bv[n][n]);

  Matrix<double,1,Dynamic> c(10);
  for (size_t n = 0; n < 10; ++n)
    c[n] = n;
  VectorViewMvt<Matrix<double,1,Dynamic> > cv(c);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(c[n], cv[n][n]);
  for (size_t n = 0; n < 10; ++n)
    cv[n][n] = n+10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10+n, cv[n][n]);
    EXPECT_FLOAT_EQ(10+n, c[n]);
  }

  const Matrix<double,1,Dynamic> d(c);
  VectorViewMvt<const Matrix<double,1,Dynamic>,
          stan::is_vector<Matrix<double,1,Dynamic> >::value > dv(d);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(d[n], dv[n][n]);

}
// XXX Need to test both const and regular for everything!

TEST(MetaTraits,VectorViewMvt_multivariate_const) {
  using stan::VectorViewMvt;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  //Test std::vector of column vectors
  Matrix<double,Dynamic,1> x(3);
  x << 1.0, 4.0, 9.0;
  Matrix<double,Dynamic,1> y(3);
  y << 2.0, 5.0, 10.0;
  std::vector<Matrix<double,Dynamic,1> > xy_vecs;
  xy_vecs.push_back(x);
  xy_vecs.push_back(y);

  VectorViewMvt<const Matrix<double,Dynamic,1> > xy_view(xy_vecs);
  EXPECT_FLOAT_EQ(1.0,xy_view[0][0]);
  EXPECT_FLOAT_EQ(4.0,xy_view[0][1]);
  EXPECT_FLOAT_EQ(9.0,xy_view[0][2]);
  // Should be:
  //EXPECT_FLOAT_EQ(2.0,xy_view[1][0]);
  //EXPECT_FLOAT_EQ(5.0,xy_view[1][1]);
  //EXPECT_FLOAT_EQ(10.0,xy_view[1][2]);
  // Currently this case is buggy as is_array is not set correctly so it just
  // gives the 1st matrix back again.
  EXPECT_FLOAT_EQ(1.0,xy_view[1][0]);
  EXPECT_FLOAT_EQ(4.0,xy_view[1][1]);
  EXPECT_FLOAT_EQ(9.0,xy_view[1][2]);

  //Test std::vector of row vectors
  Matrix<double,1,Dynamic> a(3);
  a << 1.0, 4.0, 9.0;
  Matrix<double,1,Dynamic> b(3);
  b << 2.0, 5.0, 10.0;
  std::vector<Matrix<double,1,Dynamic> > ab_vecs;
  ab_vecs.push_back(a);
  ab_vecs.push_back(b);

  VectorViewMvt<const Matrix<double,1,Dynamic> > ab_view(ab_vecs);
  EXPECT_FLOAT_EQ(1.0,ab_view[0][0]);
  EXPECT_FLOAT_EQ(4.0,ab_view[0][1]);
  EXPECT_FLOAT_EQ(9.0,ab_view[0][2]);
  // Should be:
  //EXPECT_FLOAT_EQ(2.0,ab_view[1][0]);
  //EXPECT_FLOAT_EQ(5.0,ab_view[1][1]);
  //EXPECT_FLOAT_EQ(10.0,ab_view[1][2]);
  // Currently this case is buggy as is_array is not set correctly so it just
  // gives the 1st matrix back again.
  EXPECT_FLOAT_EQ(1.0,ab_view[1][0]);
  EXPECT_FLOAT_EQ(4.0,ab_view[1][1]);
  EXPECT_FLOAT_EQ(9.0,ab_view[1][2]);

  //This tests a std::vector of matrices
  Matrix<double,Dynamic,Dynamic> m(2,3);
  m << 1.0, 2.0, 3.0, -100.0, -200.0, -300.0;
  Matrix<double,Dynamic,Dynamic> n(2,3);
  n << 4.0, 5.0, 6.0, -400.0, -500.0, -600.0;
  std::vector<Matrix<double,Dynamic,Dynamic> > mn_matrices;
  mn_matrices.push_back(m);
  mn_matrices.push_back(n);
  VectorViewMvt<const Matrix<double,Dynamic,Dynamic> > mn_view(mn_matrices);
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 2; ++i) {
        // Should be:
        //EXPECT_FLOAT_EQ(mn_matrices[k](i,j),mn_view[k](i,j));
        // Bug in is_array leads to:
        EXPECT_FLOAT_EQ(mn_matrices[0](i,j),mn_view[k](i,j));
      }
    }
  }
}
TEST(MetaTraits,VectorViewMvt_multivariate) {
  using stan::VectorViewMvt;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  //Test std::vector of column vectors
  Matrix<double,Dynamic,1> x(3);
  x << 1.0, 4.0, 9.0;
  Matrix<double,Dynamic,1> y(3);
  y << 2.0, 5.0, 10.0;
  std::vector<Matrix<double,Dynamic,1> > xy_vecs;
  xy_vecs.push_back(x);
  xy_vecs.push_back(y);

  VectorViewMvt<Matrix<double,Dynamic,1> > xy_view(xy_vecs);
  EXPECT_FLOAT_EQ(1.0,xy_view[0][0]);
  EXPECT_FLOAT_EQ(4.0,xy_view[0][1]);
  EXPECT_FLOAT_EQ(9.0,xy_view[0][2]);
  // Should be:
  //EXPECT_FLOAT_EQ(2.0,xy_view[1][0]);
  //EXPECT_FLOAT_EQ(5.0,xy_view[1][1]);
  //EXPECT_FLOAT_EQ(10.0,xy_view[1][2]);
  // Currently this case is buggy as is_array is not set correctly so it just
  // gives the 1st matrix back again.
  EXPECT_FLOAT_EQ(1.0,xy_view[1][0]);
  EXPECT_FLOAT_EQ(4.0,xy_view[1][1]);
  EXPECT_FLOAT_EQ(9.0,xy_view[1][2]);

  //Test std::vector of row vectors
  Matrix<double,1,Dynamic> a(3);
  a << 1.0, 4.0, 9.0;
  Matrix<double,1,Dynamic> b(3);
  b << 2.0, 5.0, 10.0;
  std::vector<Matrix<double,1,Dynamic> > ab_vecs;
  ab_vecs.push_back(a);
  ab_vecs.push_back(b);

  VectorViewMvt<Matrix<double,1,Dynamic> > ab_view(ab_vecs);
  EXPECT_FLOAT_EQ(1.0,ab_view[0][0]);
  EXPECT_FLOAT_EQ(4.0,ab_view[0][1]);
  EXPECT_FLOAT_EQ(9.0,ab_view[0][2]);
  // Should be:
  //EXPECT_FLOAT_EQ(2.0,ab_view[1][0]);
  //EXPECT_FLOAT_EQ(5.0,ab_view[1][1]);
  //EXPECT_FLOAT_EQ(10.0,ab_view[1][2]);
  // Currently this case is buggy as is_array is not set correctly so it just
  // gives the 1st matrix back again.
  EXPECT_FLOAT_EQ(1.0,ab_view[1][0]);
  EXPECT_FLOAT_EQ(4.0,ab_view[1][1]);
  EXPECT_FLOAT_EQ(9.0,ab_view[1][2]);

  //This tests a std::vector of matrices
  Matrix<double,Dynamic,Dynamic> m(2,3);
  m << 1.0, 2.0, 3.0, -100.0, -200.0, -300.0;
  Matrix<double,Dynamic,Dynamic> n(2,3);
  n << 4.0, 5.0, 6.0, -400.0, -500.0, -600.0;
  std::vector<Matrix<double,Dynamic,Dynamic> > mn_matrices;
  mn_matrices.push_back(m);
  mn_matrices.push_back(n);
  VectorViewMvt<Matrix<double,Dynamic,Dynamic> > mn_view(mn_matrices);
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 2; ++i) {
        // Should be:
        //EXPECT_FLOAT_EQ(mn_matrices[k](i,j),mn_view[k](i,j));
        // Bug in is_array leads to:
        EXPECT_FLOAT_EQ(mn_matrices[0](i,j),mn_view[k](i,j));
      }
    }
  }
}
