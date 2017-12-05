#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <vector>

void test_log1m_inv_logit(const double x) {
  using stan::math::var;
  using stan::math::log1m_inv_logit;
  using std::log;
  using stan::math::inv_logit;


  // test gradient
  AVEC x1 = createAVEC(x);
  AVAR f1 = log1m_inv_logit(x1[0]);
  std::vector<double> grad_f1;
  f1.grad(x1, grad_f1);

  AVEC x2 = createAVEC(x);
  AVAR f2 = log(1.0 - inv_logit(x2[0]));
  std::vector<double> grad_f2;
  f2.grad(x2, grad_f2);

  EXPECT_EQ(1U, grad_f1.size());
  EXPECT_EQ(1U, grad_f2.size());
  EXPECT_FLOAT_EQ(grad_f2[0], grad_f1[0]);

  // test value
  EXPECT_FLOAT_EQ(log(1.0 - inv_logit(x)),
                  log1m_inv_logit(var(x)).val());
}

TEST(AgradRev, log1m_inv_logit) {
  test_log1m_inv_logit(-7.2);
  test_log1m_inv_logit(0.0);
  test_log1m_inv_logit(1.9);
}

struct log1m_inv_logit_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return stan::math::log1m_inv_logit(arg1);
  }
};

TEST(AgradRev, log1m_inv_logit_NaN) {
  log1m_inv_logit_fun log1m_inv_logit_;
  test_nan(log1m_inv_logit_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  stan::math::var x = -7.2;
  test::check_varis_on_stack(stan::math::log1m_inv_logit(x));
}
