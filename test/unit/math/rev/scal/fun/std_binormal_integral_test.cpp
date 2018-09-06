#include <stan/math/prim/scal/prob/normal_cdf.hpp>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <gtest/gtest.h>
#include <limits>

using Eigen::Dynamic;
using Eigen::Matrix;

struct VVV_std_binorm_integral {
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    return stan::math::std_binormal_integral(inp_vec(0), inp_vec(1),
                                             inp_vec(2));
  }
};

template <class F1>
void compare_grad_known(const F1& f1,
                        const Matrix<double, Dynamic, 1>& exp_grad,
                        const Matrix<double, Dynamic, 1>& inp_vec) {
  double f1_eval;
  Matrix<double, Dynamic, 1> grad_f1;
  stan::math::gradient(f1, inp_vec, f1_eval, grad_f1);
  for (int i = 0; i < grad_f1.size(); ++i)
    EXPECT_FLOAT_EQ(exp_grad(i), grad_f1(i));
}

TEST(MathFunctions, binormal_integral_using) {
  using stan::math::std_binormal_integral;
}

TEST(MathFunctions, binormal_integral_throw_RV_1_nan_vvv) {
  using stan::math::var;
  var nan = std::numeric_limits<double>::quiet_NaN();
  var rho = 0.3;
  var a = nan;
  var b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan_vvv) {
  using stan::math::var;
  var nan = std::numeric_limits<var>::quiet_NaN();
  var rho = 0.3;
  var a = 2;
  var b = nan;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan_vvv) {
  using stan::math::var;
  var nan = std::numeric_limits<double>::quiet_NaN();
  var rho = nan;
  var a = 2;
  var b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_neg_vvv) {
  using stan::math::var;
  var rho = -1.3;
  var a = 2;
  var b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_gt_one_vvv) {
  using stan::math::var;
  var rho = 1.3;
  var a = 2;
  var b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_no_throw_vvv) {
  using stan::math::var;
  var rho = 0.3;
  var a = 2;
  var b = 1;
  EXPECT_NO_THROW(stan::math::std_binormal_integral(a, b, rho));
}
TEST(MathFunctions, binormal_integral_val_boundaries_test) {
  // Independent normal RVs
  using stan::math::var;
  var rho = 0;
  var a = -0.4;
  var b = 2.7;
  var f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(a.val()) * stan::math::Phi(b.val()), f.val());

  // Perfectly correlated RVs
  rho = 1;
  a = -3.4;
  b = 3.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(a.val()), f.val());

  // Perfectly anticorrelated RVs
  rho = -1;
  a = 2.4;
  b = 1.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(a.val()) + stan::math::Phi(b.val()) - 1,
                  f.val());

  // Perfectly anticorrelated RVs
  rho = -1;
  a = -2.4;
  b = 1.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0, f.val());

  // a = rho * b
  rho = -0.7;
  b = 1.7;
  a = rho * b;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.5 / stan::math::pi() * std::exp(-0.5 * b.val() * b.val())
                          * std::asin(rho.val())
                      + stan::math::Phi(a.val()) * stan::math::Phi(b.val()),
                  f.val());
  // b = rho * a
  rho = -0.7;
  a = 1.7;
  b = rho * a;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.5 / stan::math::pi() * std::exp(-0.5 * a.val() * a.val())
                          * std::asin(rho.val())
                      + stan::math::Phi(a.val()) * stan::math::Phi(b.val()),
                  f.val());
  rho = 0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(1, f.val());

  rho = 0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0, f.val());

  rho = -0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0, f.val());

  rho = -0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(1, f.val());

  rho = -0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(1.5), f.val());

  rho = 0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(1.5), f.val());

  rho = 0.7;
  b = 2.5;
  a = std::numeric_limits<double>::infinity();
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(2.5), f.val());

  rho = -0.7;
  b = 0.5;
  a = std::numeric_limits<double>::infinity();
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(0.5), f.val());
}
TEST(MathFunctions, binormal_integral_val_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,
  // b), corr = matrix(c(1, rho, rho, 1), 2, 2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  using stan::math::var;
  var rho = 0.3;
  var a = -0.4;
  var b = 2.7;
  var f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.344276561500873, f.val());

  rho = 0.99;
  a = -0.4;
  b = 2.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.3445782583896758, f.val());

  rho = 0.99;
  a = 2.5;
  b = 2.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.9937227710497979, f.val());

  rho = 0.99;
  a = 3.5;
  b = 3.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.9997643606337163, f.val());

  rho = -0.99;
  a = -4.5;
  b = 4.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(2.146032113348184e-06, f.val());

  rho = -0.99;
  a = -4.5;
  b = 10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(3.397673124738709e-06, f.val());

  rho = 0.99;
  a = 4.5;
  b = -6.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(4.016000583859118e-11, f.val());

  rho = 0.5;
  a = -4.5;
  b = -10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(5.612932952882046e-24, f.val());

  rho = 0.99;
  a = -4.5;
  b = -10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(7.619853024160583e-24, f.val());

  rho = 0.5;
  a = -4.5;
  b = -10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(5.612932952882069e-24, f.val());

  rho = -0.3;
  a = -3.3;
  b = -3.4;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(-21.05454315638917, log(f.val()));
}
TEST(MathFunctions, binormal_integral_grad_test_vvv_owens) {
  using stan::math::normal_cdf;
  using stan::math::var;
  using std::exp;
  using std::sqrt;
  var rho = 0.3;
  var a = -0.4;
  var b = 2.7;
  var f = stan::math::std_binormal_integral(a, b, rho);

  double gf_1 = 1 / sqrt(2 * stan::math::pi()) * exp(-0.5 * a.val() * a.val())
                * normal_cdf(b.val(), rho.val() * a.val(),
                             sqrt(1 - rho.val() * rho.val()));
  double gf_2 = 1 / sqrt(2 * stan::math::pi()) * exp(-0.5 * b.val() * b.val())
                * normal_cdf(a.val(), rho.val() * b.val(),
                             sqrt(1 - rho.val() * rho.val()));
  double gf_3 = 0.5 / (stan::math::pi() * sqrt(1 - rho.val() * rho.val()))
                * exp(-0.5 / (1 - rho.val() * rho.val())
                          * (a.val() - rho.val() * b.val())
                          * (a.val() - rho.val() * b.val())
                      - 0.5 * b.val() * b.val());
  AVEC x = createAVEC(a, b, rho);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(gf_1, grad_f[0]);
  EXPECT_FLOAT_EQ(gf_2, grad_f[1]);
  EXPECT_FLOAT_EQ(gf_3, grad_f[2]);
}
TEST(MathFunctions, binormal_integral_grad_test_vvv_tanh_sinh) {
  using stan::math::normal_cdf;
  using stan::math::var;
  using std::exp;
  using std::sqrt;
  var rho = -0.5;
  var a = -2.3;
  var b = -4.4;
  var f = stan::math::std_binormal_integral(a, b, rho);

  double gf_1 = 1 / sqrt(2 * stan::math::pi()) * exp(-0.5 * a.val() * a.val())
                * normal_cdf(b.val(), rho.val() * a.val(),
                             sqrt(1 - rho.val() * rho.val()));
  double gf_2 = 1 / sqrt(2 * stan::math::pi()) * exp(-0.5 * b.val() * b.val())
                * normal_cdf(a.val(), rho.val() * b.val(),
                             sqrt(1 - rho.val() * rho.val()));
  double gf_3 = 0.5 / (stan::math::pi() * sqrt(1 - rho.val() * rho.val()))
                * exp(-0.5 / (1 - rho.val() * rho.val())
                          * (a.val() - rho.val() * b.val())
                          * (a.val() - rho.val() * b.val())
                      - 0.5 * b.val() * b.val());
  AVEC x = createAVEC(a, b, rho);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(gf_1, grad_f[0]);
  EXPECT_FLOAT_EQ(gf_2, grad_f[1]);
  EXPECT_FLOAT_EQ(gf_3, grad_f[2]);
}
TEST(MathFunctions, binormal_integral_grad_test_vvv_at_boundary) {
  VVV_std_binorm_integral test_fun;
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << -std::numeric_limits<double>::infinity(), 2, 0.5;
  Matrix<double, Dynamic, 1> expected_grad(3);
  expected_grad << 0, 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << 2, -std::numeric_limits<double>::infinity(), 0.5;
  expected_grad << 0, 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << 2, std::numeric_limits<double>::infinity(), 0.5;
  expected_grad << exp(stan::math::std_normal_lpdf(2)), 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << std::numeric_limits<double>::infinity(), 1, 0.5;
  expected_grad << 0, exp(stan::math::std_normal_lpdf(1)), 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << std::numeric_limits<double>::infinity(), 1, 1;
  expected_grad << 0, exp(stan::math::std_normal_lpdf(1)), 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << std::numeric_limits<double>::infinity(), 1, -1;
  expected_grad << 0, exp(stan::math::std_normal_lpdf(1)), 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), 0.4;
  expected_grad << 0, 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(), 0.4;
  expected_grad << 0, 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << -std::numeric_limits<double>::infinity(), 1, 1;
  expected_grad << 0, 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << 1, -std::numeric_limits<double>::infinity(), 1;
  expected_grad << 0, 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << 1, -std::numeric_limits<double>::infinity(), -1;
  expected_grad << 0, 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << 2, 3, 1;
  expected_grad << 1 / sqrt(2 * stan::math::pi()) * exp(-0.5 * 4), 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << 3, 1, 1;
  expected_grad << 0, 1 / sqrt(2 * stan::math::pi()) * exp(-0.5), 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << 3, 1, -1;
  expected_grad << 1 / sqrt(2 * stan::math::pi()) * exp(-0.5 * 9),
      1 / sqrt(2 * stan::math::pi()) * exp(-0.5), 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);

  inp_vec << -1, 1, -1;
  expected_grad << 0, 0, 0;
  compare_grad_known(test_fun, expected_grad, inp_vec);
}
TEST(MathFunctions, binormal_integral_grad_test_vvv_tan_sinh_near_boundary) {
  using stan::math::normal_cdf;
  using stan::math::var;
  using std::exp;
  using std::sqrt;
  var rho = 0.99;
  var a = -5.3;
  var b = -5.4;
  var f = stan::math::std_binormal_integral(a, b, rho);

  double one_minus_r_sq = (1 + rho.val()) * (1 - rho.val());
  double gf_1
      = 1 / sqrt(2 * stan::math::pi()) * exp(-0.5 * a.val() * a.val())
        * normal_cdf(b.val(), rho.val() * a.val(), sqrt(one_minus_r_sq));
  double gf_2
      = 1 / sqrt(2 * stan::math::pi()) * exp(-0.5 * b.val() * b.val())
        * normal_cdf(a.val(), rho.val() * b.val(), sqrt(one_minus_r_sq));
  double gf_3 = 0.5 / (stan::math::pi() * sqrt(one_minus_r_sq))
                * exp(-0.5 / one_minus_r_sq * (a.val() - rho.val() * b.val())
                          * (a.val() - rho.val() * b.val())
                      - 0.5 * b.val() * b.val());
  AVEC x = createAVEC(a, b, rho);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(gf_1, grad_f[0]);
  EXPECT_FLOAT_EQ(gf_2, grad_f[1]);
  EXPECT_FLOAT_EQ(gf_3, grad_f[2]);
}

TEST(AgradRev, check_varis_on_stack) {
  stan::math::var rho = 0.5;
  stan::math::var a = 2.0;
  stan::math::var b = 1.0;
  test::check_varis_on_stack(stan::math::std_binormal_integral(a, b, rho));
  test::check_varis_on_stack(
      stan::math::std_binormal_integral(a.val(), b, rho));
  test::check_varis_on_stack(
      stan::math::std_binormal_integral(a, b.val(), rho));
  test::check_varis_on_stack(
      stan::math::std_binormal_integral(a.val(), b.val(), rho));
  test::check_varis_on_stack(
      stan::math::std_binormal_integral(a, b, rho.val()));
  test::check_varis_on_stack(
      stan::math::std_binormal_integral(a.val(), b, rho.val()));
  test::check_varis_on_stack(
      stan::math::std_binormal_integral(a, b.val(), rho.val()));
}
