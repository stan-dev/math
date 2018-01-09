#include <stan/math/rev/scal.hpp>
#include <stan/math/prim/scal/fun/grad_reg_lower_inc_gamma.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/gamma.hpp>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev, gamma_p_var_var) {
  AVAR a = 0.5001;
  AVAR b = 1.0001;
  AVAR f = gamma_p(a, b);
  EXPECT_FLOAT_EQ(boost::math::gamma_p(0.5001, 1.0001), f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(-0.3898178624664172, g[0]);
  EXPECT_FLOAT_EQ(boost::math::gamma_p_derivative(0.5001, 1.0001), g[1]);

  a = -0.5001;
  EXPECT_THROW(gamma_p(a, b), std::domain_error);

  b = -1.0;
  EXPECT_THROW(gamma_p(a, b), std::domain_error);
}
TEST(AgradRev, gamma_p_double_var) {
  double a = 0.5001;
  AVAR b = 1.0001;
  AVAR f = gamma_p(a, b);
  EXPECT_FLOAT_EQ(boost::math::gamma_p(0.5001, 1.0001), f.val());

  AVEC x = createAVEC(b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(boost::math::gamma_p_derivative(0.5001, 1.0001), g[0]);

  a = -0.5;
  EXPECT_THROW(gamma_p(a, b), std::domain_error);

  b = -1.0;
  EXPECT_THROW(gamma_p(a, b), std::domain_error);
}
TEST(AgradRev, gamma_p_var_double) {
  AVAR a = 0.5001;
  double b = 1.0001;
  AVAR f = gamma_p(a, b);
  EXPECT_FLOAT_EQ(boost::math::gamma_p(0.5001, 1.0001), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(-0.3898178624664172, g[0]);

  a = -0.5;
  EXPECT_THROW(gamma_p(a, b), std::domain_error);

  b = -1.0;
  EXPECT_THROW(gamma_p(a, b), std::domain_error);
}

struct gamma_p_fun {
  template <typename T0, typename T1>
  inline typename stan::return_type<T0, T1>::type operator()(
      const T0& arg1, const T1& arg2) const {
    return gamma_p(arg1, arg2);
  }
};

TEST(AgradRev, gamma_p_nan) {
  gamma_p_fun gamma_p_;
  test_nan(gamma_p_, 0.5001, 1.0001, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 0.5001;
  AVAR b = 1.0001;
  test::check_varis_on_stack(stan::math::gamma_p(a, b));
  test::check_varis_on_stack(stan::math::gamma_p(a, 1.0001));
  test::check_varis_on_stack(stan::math::gamma_p(0.5001, b));
}
