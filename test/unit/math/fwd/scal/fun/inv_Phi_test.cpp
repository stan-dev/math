#include <gtest/gtest.h>
#include <stan/math/fwd/scal/fun/Phi.hpp>
#include <stan/math/fwd/scal/fun/inv_Phi.hpp>
#include <stan/math/prim/scal/fun/inv_Phi.hpp>
#include <stan/math/prim/scal/prob/normal_log.hpp>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal/fun/exp.hpp>

TEST(AgradFwdinv_Phi,Fvar) {
  using stan::math::fvar;
  using stan::math::inv_Phi;
  using stan::math::Phi;

  fvar<double> x = 0.1;
  x.d_ = 1.0;

  fvar<double> y = Phi(inv_Phi(x));

  EXPECT_FLOAT_EQ(0.1, y.val_);
  EXPECT_FLOAT_EQ(1.0, y.d_);
}
TEST(AgradFwdinv_Phi, FvarFvarDouble) {
  using stan::math::fvar;
  using stan::math::inv_Phi;
  using stan::math::Phi;

  fvar<fvar<double> > x;
  x.val_.val_ = 0.1;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > a = Phi(inv_Phi(x));

  EXPECT_FLOAT_EQ(0.1, a.val_.val_);
  EXPECT_FLOAT_EQ(1.0, a.val_.d_);
  EXPECT_FLOAT_EQ(0, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);
}

struct inv_Phi_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return inv_Phi(arg1);
  }
};

TEST(AgradFwdinv_Phi,inv_Phi_NaN) {
  inv_Phi_fun inv_Phi_;
  test_nan_fwd(inv_Phi_,true);
}
