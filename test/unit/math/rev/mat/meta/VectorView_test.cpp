#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorView_matrix_var) {
  using stan::VectorView;
  using stan::math::var;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  
  Matrix<var,Dynamic,1> a(10);
  for (size_t n = 0; n < 10; ++n) 
    a[n] = n;
  VectorView<Matrix<var,Dynamic,1> > av(a);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(a[n].val(), av[n].val());
  for (size_t n = 0; n < 10; ++n)
    av[n] = n+10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10+n, av[n].val());
    EXPECT_FLOAT_EQ(10+n, a[n].val());
  }

  const Matrix<var,Dynamic,1> b(a);
  VectorView<const Matrix<var,Dynamic,1> > bv(b);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(b[n].val(), bv[n].val());
  
  Matrix<var,1,Dynamic> c(10);
  for (size_t n = 0; n < 10; ++n) 
    c[n] = n;
  VectorView<Matrix<var,1,Dynamic> > cv(c);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(c[n].val(), cv[n].val());
  for (size_t n = 0; n < 10; ++n) 
    cv[n] = n+10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(10+n, cv[n].val());
    EXPECT_FLOAT_EQ(10+n, c[n].val());
  }

  const Matrix<var,1,Dynamic> d(c);
  VectorView<const Matrix<var,1,Dynamic> > dv(d);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(d[n].val(), dv[n].val());
}
