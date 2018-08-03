#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/rev/arr/fun/util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev, log_sum_exp_vv) {
  AVAR a = 5.0;
  AVAR b = 2.0;
  AVAR f = log_sum_exp(a, b);
  EXPECT_FLOAT_EQ(std::log(std::exp(5) + std::exp(2)), f.val());

  AVEC x = createAVEC(a, b);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(std::exp(5.0) / (std::exp(5.0) + std::exp(2.0)), grad_f[0]);
  EXPECT_FLOAT_EQ(std::exp(2.0) / (std::exp(5.0) + std::exp(2.0)), grad_f[1]);

  // underflow example
  a = 1000;
  b = 10;
  f = log_sum_exp(a, b);
  EXPECT_FLOAT_EQ(std::log(std::exp(0.0) + std::exp(-990.0)) + 1000.0, f.val());

  x = createAVEC(a, b);
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(
      std::exp(1000.0 - (std::log(std::exp(0.0) + std::exp(-999.0)) + 1000)),
      grad_f[0]);
  EXPECT_FLOAT_EQ(
      std::exp(10.0 - (std::log(std::exp(0.0) + std::exp(-999.0)) + 1000)),
      grad_f[1]);
}
TEST(AgradRev, log_sum_exp_vd) {
  AVAR a = 5.0;
  double b = 2.0;
  AVAR f = log_sum_exp(a, b);
  EXPECT_FLOAT_EQ(std::log(std::exp(5) + std::exp(2)), f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(std::exp(5.0) / (std::exp(5.0) + std::exp(2.0)), grad_f[0]);

  // underflow example
  a = 1000;
  b = 10;
  f = log_sum_exp(a, b);
  EXPECT_FLOAT_EQ(std::log(std::exp(0.0) + std::exp(-990.0)) + 1000.0, f.val());

  x = createAVEC(a);
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(
      std::exp(1000.0 - (std::log(std::exp(0.0) + std::exp(-999.0)) + 1000)),
      grad_f[0]);
}
TEST(AgradRev, log_sum_exp_dv) {
  double a = 5.0;
  AVAR b = 2.0;
  AVAR f = log_sum_exp(a, b);
  EXPECT_FLOAT_EQ(std::log(std::exp(5) + std::exp(2)), f.val());

  AVEC x = createAVEC(b);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(std::exp(2.0) / (std::exp(5.0) + std::exp(2.0)), grad_f[0]);

  // underflow example
  a = 10;
  b = 1000;
  f = log_sum_exp(a, b);
  EXPECT_FLOAT_EQ(std::log(std::exp(0.0) + std::exp(-990.0)) + 1000.0, f.val());

  x = createAVEC(b);
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(
      std::exp(1000.0 - (std::log(std::exp(0.0) + std::exp(-999.0)) + 1000)),
      grad_f[0]);
}

void test_log_sum_exp_2_vv(double a_val, double b_val) {
  using stan::math::log_sum_exp;
  using std::exp;
  using std::log;

  AVAR a(a_val);
  AVAR b(b_val);

  AVEC x = createAVEC(a, b);
  AVAR f = log_sum_exp(a, b);
  VEC g;
  f.grad(x, g);

  double f_val = f.val();

  stan::math::var a2(a_val);
  stan::math::var b2(b_val);
  AVEC x2 = createAVEC(a2, b2);
  AVAR f2 = log(exp(a2) + exp(b2));
  VEC g2;
  f2.grad(x2, g2);

  EXPECT_FLOAT_EQ(f2.val(), f_val);
  EXPECT_EQ(2U, g.size());
  EXPECT_EQ(2U, g2.size());
  EXPECT_FLOAT_EQ(g2[0], g[0]);
  EXPECT_FLOAT_EQ(g2[1], g[1]);
}
void test_log_sum_exp_2_vd(double a_val, double b) {
  using stan::math::log_sum_exp;
  using std::exp;
  using std::log;

  AVAR a(a_val);
  AVEC x = createAVEC(a);
  AVAR f = log_sum_exp(a, b);
  VEC g;
  f.grad(x, g);

  double f_val = f.val();

  stan::math::var a2(a_val);
  AVEC x2 = createAVEC(a2);
  AVAR f2 = log(exp(a2) + exp(b));
  VEC g2;
  f2.grad(x2, g2);

  EXPECT_FLOAT_EQ(f2.val(), f_val);
  EXPECT_EQ(1U, g.size());
  EXPECT_EQ(1U, g2.size());
  EXPECT_FLOAT_EQ(g2[0], g[0]);
}
void test_log_sum_exp_2_dv(double a, double b_val) {
  using stan::math::log_sum_exp;
  using std::exp;
  using std::log;

  AVAR b(b_val);
  AVEC x = createAVEC(b);
  AVAR f = log_sum_exp(a, b);
  VEC g;
  f.grad(x, g);

  double f_val = f.val();

  AVAR b2(b_val);
  AVEC x2 = createAVEC(b2);
  AVAR f2 = log(exp(a) + exp(b2));
  VEC g2;
  f2.grad(x2, g2);

  EXPECT_FLOAT_EQ(f2.val(), f_val);
  EXPECT_EQ(1U, g.size());
  EXPECT_EQ(1U, g2.size());
  EXPECT_FLOAT_EQ(g2[0], g[0]);
}

void test_log_sum_exp_2(double a, double b) {
  test_log_sum_exp_2_vv(a, b);
  test_log_sum_exp_2_vd(a, b);
  test_log_sum_exp_2_dv(a, b);
}

TEST(AgradRev, log_sum_exp_2) {
  test_log_sum_exp_2(0.0, 0.0);
  test_log_sum_exp_2(1.0, 2.0);
  test_log_sum_exp_2(2.0, 1.0);
  test_log_sum_exp_2(-2.0, 15.0);
  test_log_sum_exp_2(2.0, -15.0);
}

struct log_sum_exp_fun {
  template <typename T0, typename T1>
  inline typename stan::return_type<T0, T1>::type operator()(
      const T0& arg1, const T1& arg2) const {
    return log_sum_exp(arg1, arg2);
  }
};

TEST(AgradRev, log_sum_exp_nan) {
  log_sum_exp_fun log_sum_exp_;
  test_nan(log_sum_exp_, 3.0, 5.0, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 5.0;
  AVAR b = 2.0;
  test::check_varis_on_stack(log_sum_exp(a, b));
  test::check_varis_on_stack(log_sum_exp(5.0, b));
  test::check_varis_on_stack(log_sum_exp(a, 2.0));
}
