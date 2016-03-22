#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/type_traits.hpp>

TEST(MetaTraits, VectorView_double)  {
  using stan::VectorView;

  double d(10);
  VectorView<double> dv(d);
  EXPECT_FLOAT_EQ(d, dv[0]);
  dv[1] = 2.0;
  EXPECT_FLOAT_EQ(2.0, dv[0]);
  EXPECT_FLOAT_EQ(2.0, d);

  const double c(10);
  VectorView<const double> cv(c);
  EXPECT_FLOAT_EQ(c, cv[0]);
}

TEST(MetaTraits,VectorView) {
  using stan::VectorView;

  double x = 5;
  VectorView<const double> x_VectorView(x);
  EXPECT_FLOAT_EQ(x,x_VectorView[0]);
  EXPECT_FLOAT_EQ(x,x_VectorView[1]);
  EXPECT_FLOAT_EQ(x,x_VectorView[2]);
}

TEST(MetaTraits, VectorView_double_star) {
  using stan::VectorView;
  double a[10];
  double *a_star = &a[0];
  for (size_t n = 0; n < 10; ++n)
    a[n] = n;
  VectorView<double*,true> av(a_star);
  for (size_t n = 0; n < 10; ++n) 
    EXPECT_FLOAT_EQ(a[n], av[n]);
  for (size_t n = 0; n < 10; ++n)
    av[n] = n+10;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(n+10, a[n]);
    EXPECT_FLOAT_EQ(n+10, av[n]);
  }

  double b(20);
  double *b_star = &b;
  VectorView<double*,false> bv(b_star);
  for (size_t n = 0; n < 10; ++n) 
    EXPECT_FLOAT_EQ(20, bv[n]);
  bv[1] = 10;
  EXPECT_FLOAT_EQ(10, bv[0]);
  EXPECT_FLOAT_EQ(10, b);
}

TEST(MetaTraits, VectorView_throw_if_accessed) {
  using stan::VectorView;
  double d(10);
  VectorView<double, false, true> dv(d);
  EXPECT_THROW(dv[0], std::logic_error);
}
