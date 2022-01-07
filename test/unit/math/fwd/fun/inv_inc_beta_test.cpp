#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdMatrixIncBetaInv, fd_scalar) {
  using stan::math::fvar;
  using stan::math::inv_inc_beta;
  fvar<double> a = 6;
  fvar<double> b = 2;
  fvar<double> p = 0.9;
  a.d_ = 1.0;
  b.d_ = 1.0;
  p.d_ = 1.0;

  fvar<double> res = inv_inc_beta(a, b, p);

  EXPECT_FLOAT_EQ(res.d_, 0.0117172527399 - 0.0680999818473 + 0.455387298585);
}

TEST(AgradFwdMatrixIncBetaInv, ffd_scalar) {
  using stan::math::fvar;
  using stan::math::inv_inc_beta;
  fvar<fvar<double>> a = 7;
  fvar<fvar<double>> b = 4;
  fvar<fvar<double>> p = 0.15;
  a.val_.d_ = 1.0;
  b.val_.d_ = 1.0;
  p.val_.d_ = 1.0;

  fvar<fvar<double>> res = inv_inc_beta(a, b, p);

  EXPECT_FLOAT_EQ(res.val_.d_,
                  0.0428905418857 - 0.0563420377808 + 0.664919819507);
}
