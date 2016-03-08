#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
TEST(MetaTraits, VectorView_matrix_double) {
  using stan::VectorView;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  
  Matrix<double,Dynamic,1> a(10);
  for (size_t n = 0; n < 10; ++n) 
    a[n] = n;
  VectorView<Matrix<double,Dynamic,1> > av(a);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(a[n], av[n]);
  for (size_t n = 0; n < 10; ++n)
    av[n] = n+10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10+n, av[n]);
    EXPECT_FLOAT_EQ(10+n, a[n]);
  }

  const Matrix<double,Dynamic,1> b(a);
  VectorView<const Matrix<double,Dynamic,1> > bv(b);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(b[n], bv[n]);
  
  Matrix<double,1,Dynamic> c(10);
  for (size_t n = 0; n < 10; ++n) 
    c[n] = n;
  VectorView<Matrix<double,1,Dynamic> > cv(c);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(c[n], cv[n]);
  for (size_t n = 0; n < 10; ++n) 
    cv[n] = n+10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10+n, cv[n]);
    EXPECT_FLOAT_EQ(10+n, c[n]);
  }

  const Matrix<double,1,Dynamic> d(c);
  VectorView<const Matrix<double,1,Dynamic>, 
          stan::is_vector<Matrix<double,1,Dynamic> >::value > dv(d);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(d[n], dv[n]);
}

TEST(MetaTraits,VectorView) {
  using stan::VectorView;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<double,Dynamic,1> v(3);
  v << 1.0, 4.0, 9.0;
  VectorView<const Matrix<double,Dynamic,1> > v_VectorView(v);
  EXPECT_FLOAT_EQ(1.0,v_VectorView[0]);
  EXPECT_FLOAT_EQ(4.0,v_VectorView[1]);
  EXPECT_FLOAT_EQ(9.0,v_VectorView[2]);

  Matrix<double,1,Dynamic> rv(3);
  rv << 1.0, 4.0, 9.0;
  VectorView<const Matrix<double,1,Dynamic> > rv_VectorView(rv);
  EXPECT_FLOAT_EQ(1.0,rv_VectorView[0]);
  EXPECT_FLOAT_EQ(4.0,rv_VectorView[1]);
  EXPECT_FLOAT_EQ(9.0,rv_VectorView[2]);

  Matrix<double,Dynamic,Dynamic> m(2,3);
  m << 1.0, 2.0, 3.0, -100.0, -200.0, -300.0;
  VectorView<const Matrix<double,Dynamic,Dynamic> > m_VectorView(m);
  int pos = 0;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 2; ++i) {
      EXPECT_FLOAT_EQ(m(i,j),m_VectorView[pos]);
      ++pos;
    }
  }
}
