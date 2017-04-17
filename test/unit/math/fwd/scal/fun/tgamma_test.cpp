#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdTgamma, Fvar) {
  using stan::math::fvar;
  using boost::math::tgamma;
  using boost::math::digamma;

  fvar<double> x(0.5,1.0);
  fvar<double> a = tgamma(x);
  EXPECT_FLOAT_EQ(tgamma(0.5), a.val_);
  EXPECT_FLOAT_EQ(tgamma(0.5) * digamma(0.5), a.d_);
}

struct tgamma_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return tgamma(arg1);
  }
};

TEST(AgradFwdTgamma,tgamma_NaN) {
  tgamma_fun tgamma_;
  test_nan_fwd(tgamma_,false);
}
