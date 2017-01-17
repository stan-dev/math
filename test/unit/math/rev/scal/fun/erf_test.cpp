#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <boost/math/special_functions/erf.hpp>

TEST(AgradRev,erf) {
  AVAR a = 1.3;
  AVAR f = erf(a);
  EXPECT_FLOAT_EQ(boost::math::erf(1.3), f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x,grad_f);
  EXPECT_FLOAT_EQ(2.0 / std::sqrt(boost::math::constants::pi<double>()) * std::exp(- 1.3 * 1.3), grad_f[0]);
}
struct erf_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return erf(arg1);
  }
};

TEST(AgradRev,erf_NaN) {
  erf_fun erf_;
  test_nan(erf_,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 1.3;
  test::check_varis_on_stack(stan::math::erf(a));
}
