#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdMatrixIncBetaInv, fd_scalar) {
  using stan::math::fvar;
  using stan::math::inc_beta_inv;
  fvar<double> a = 1;
  fvar<double> b = 2;
  fvar<double> p = 0.5;
  a.d_ = 1.0;
  b.d_ = 1.0;
  p.d_ = 1.0;

  fvar<double> res = inc_beta_inv(a, b, p);

  EXPECT_FLOAT_EQ(res.d_, 0.287698278597 - 0.122532267934 + 0.707106781187);
}

TEST(AgradFwdMatrixIncBetaInv, ffd_scalar) {
  using stan::math::fvar;
  using stan::math::inc_beta_inv;
  fvar<fvar<double>> a = 1;
  fvar<fvar<double>> b = 2;
  fvar<fvar<double>> p = 0.5;
  a.val_.d_ = 1.0;
  b.val_.d_ = 1.0;
  p.val_.d_ = 1.0;

  fvar<fvar<double>> res = inc_beta_inv(a, b, p);

  EXPECT_FLOAT_EQ(res.val_.d_, 0.287698278597 - 0.122532267934 + 0.707106781187);
}