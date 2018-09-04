#include <test/unit/math/prim/mat/prob/binormal_copula_cdf_test_helper.hpp>
#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::var;
using std::vector;

template <class F1, class F2>
void compare_grad(const F1& f1, const F2& f2,
                  const Matrix<double, Dynamic, 1>& inp_vec) {
  double f1_eval;
  Matrix<double, Dynamic, 1> grad_f1;
  double f2_eval;
  Matrix<double, Dynamic, 1> grad_f2;
  stan::math::gradient(f1, inp_vec, f1_eval, grad_f1);
  stan::math::gradient(f2, inp_vec, f2_eval, grad_f2);
  EXPECT_FLOAT_EQ(f1_eval, f2_eval);
  for (int i = 0; i < grad_f1.size(); ++i)
    EXPECT_FLOAT_EQ(grad_f1(i), grad_f2(i)) << i;
}

template <class F1>
void compare_grad_known(const F1& f1,
                        const Matrix<double, Dynamic, 1>& expected_grad,
                        const Matrix<double, Dynamic, 1>& inp_vec) {
  double f1_eval;
  Matrix<double, Dynamic, 1> grad_f1;
  stan::math::gradient(f1, inp_vec, f1_eval, grad_f1);
  for (int i = 0; i < grad_f1.size(); ++i)
    EXPECT_FLOAT_EQ(grad_f1(i), expected_grad(i)) << i;
}

template <class F1>
void compare_fd(const F1& f1, const Matrix<double, Dynamic, 1>& inp_vec) {
  double f1_eval;
  Matrix<double, Dynamic, 1> grad_f1;
  double f1_fd;
  Matrix<double, Dynamic, 1> grad_f1_fd;
  stan::math::gradient(f1, inp_vec, f1_eval, grad_f1);
  stan::math::finite_diff_gradient(f1, inp_vec, f1_fd, grad_f1_fd);
  EXPECT_FLOAT_EQ(f1_eval, f1_fd);
  for (int i = 0; i < grad_f1.size(); ++i)
    EXPECT_FLOAT_EQ(grad_f1(i), grad_f1_fd(i)) << i;
}

void test_invalid_args(const vector<Matrix<double, Dynamic, 1>>& invalid_args) {
  for (size_t i = 0; i < invalid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = invalid_args[i];
    auto args_V_vv = make_args<var, var, Dynamic, 1>(inp);
    auto args_V_dv = make_args<double, var, Dynamic, 1>(inp);
    auto args_V_vd = make_args<var, double, Dynamic, 1>(inp);
    auto args_RV_vv = make_args<var, var, 1, Dynamic>(inp);
    auto args_RV_dv = make_args<double, var, 1, Dynamic>(inp);
    auto args_RV_vd = make_args<var, double, 1, Dynamic>(inp);
    EXPECT_THROW(call_args(args_V_vv), std::domain_error);
    EXPECT_THROW(call_args(args_V_dv), std::domain_error);
    EXPECT_THROW(call_args(args_V_vd), std::domain_error);
    EXPECT_THROW(call_args(args_RV_vv), std::domain_error);
    EXPECT_THROW(call_args(args_RV_dv), std::domain_error);
    EXPECT_THROW(call_args(args_RV_vd), std::domain_error);
  }
}

void test_bad_sizes(const vector<Matrix<double, Dynamic, 1>>& invalid_args) {
  for (size_t i = 0; i < invalid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = invalid_args[i];
    auto args_V_vv = make_args<var, var, Dynamic, 1>(inp);
    auto args_V_dv = make_args<double, var, Dynamic, 1>(inp);
    auto args_V_vd = make_args<var, double, Dynamic, 1>(inp);
    auto args_RV_vv = make_args<var, var, 1, Dynamic>(inp);
    auto args_RV_dv = make_args<double, var, 1, Dynamic>(inp);
    auto args_RV_vd = make_args<var, double, 1, Dynamic>(inp);
    EXPECT_THROW(call_args(args_V_vv), std::invalid_argument);
    EXPECT_THROW(call_args(args_V_dv), std::invalid_argument);
    EXPECT_THROW(call_args(args_V_vd), std::invalid_argument);
    EXPECT_THROW(call_args(args_RV_vv), std::invalid_argument);
    EXPECT_THROW(call_args(args_RV_dv), std::invalid_argument);
    EXPECT_THROW(call_args(args_RV_vd), std::invalid_argument);
  }
}

void test_valid_args(const vector<Matrix<double, Dynamic, 1>>& valid_args) {
  for (size_t i = 0; i < valid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = valid_args[i];
    auto args_V_vv = make_args<var, var, Dynamic, 1>(inp);
    auto args_V_dv = make_args<double, var, Dynamic, 1>(inp);
    auto args_V_vd = make_args<var, double, Dynamic, 1>(inp);
    auto args_RV_vv = make_args<var, var, 1, Dynamic>(inp);
    auto args_RV_dv = make_args<double, var, 1, Dynamic>(inp);
    auto args_RV_vd = make_args<var, double, 1, Dynamic>(inp);
    EXPECT_NO_THROW(call_args(args_V_vv));
    EXPECT_NO_THROW(call_args(args_V_dv));
    EXPECT_NO_THROW(call_args(args_V_vd));
    EXPECT_NO_THROW(call_args(args_RV_vv));
    EXPECT_NO_THROW(call_args(args_RV_dv));
    EXPECT_NO_THROW(call_args(args_RV_vd));
  }
}

void test_vals(const vector<Matrix<double, Dynamic, 1>>& valid_args,
               const vector<double>& vals) {
  for (size_t i = 0; i < valid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = valid_args[i];
    auto args_V_vv = make_args<var, var, Dynamic, 1>(inp);
    auto args_V_dv = make_args<double, var, Dynamic, 1>(inp);
    auto args_V_vd = make_args<var, double, Dynamic, 1>(inp);
    auto args_RV_vv = make_args<var, var, 1, Dynamic>(inp);
    auto args_RV_dv = make_args<double, var, 1, Dynamic>(inp);
    auto args_RV_vd = make_args<var, double, 1, Dynamic>(inp);
    EXPECT_FLOAT_EQ(value_of_rec(call_args(args_V_vv)), vals[i]);
    EXPECT_FLOAT_EQ(value_of_rec(call_args(args_V_dv)), vals[i]);
    EXPECT_FLOAT_EQ(value_of_rec(call_args(args_V_vd)), vals[i]);
    EXPECT_FLOAT_EQ(value_of_rec(call_args(args_RV_vv)), vals[i]);
    EXPECT_FLOAT_EQ(value_of_rec(call_args(args_RV_dv)), vals[i]);
    EXPECT_FLOAT_EQ(value_of_rec(call_args(args_RV_vd)), vals[i]);
  }
}

TEST(MathFunctions, binormal_copula_cdf_using) {
  using stan::math::binormal_copula_cdf;
}

TEST(MathFunctions, binormal_copula_cdf_throw_invalid_args_all_types) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  Matrix<double, Dynamic, 1> inp(3);
  vector<Matrix<double, Dynamic, 1>> invalid_args;
  // RV1 nan
  inp << nan, 0.3, 0.3;
  invalid_args.push_back(inp);
  // RV2 nan
  inp << 0.3, nan, 0.3;
  invalid_args.push_back(inp);
  // rho nan
  inp << 0.3, 0.3, nan;
  invalid_args.push_back(inp);
  // RV1 gt 1
  inp << 1.2, 0.3, 0.4;
  invalid_args.push_back(inp);
  // RV1 lt 0
  inp << -1.2, 0.3, 0.4;
  invalid_args.push_back(inp);
  // RV2 gt 1
  inp << 0.3, 1.2, 0.4;
  invalid_args.push_back(inp);
  // RV2 lt 0
  inp << 0.2, -1.2, 0.3;
  invalid_args.push_back(inp);
  // rho lt neg 1
  inp << 0.3, 0.3, -1.3;
  invalid_args.push_back(inp);
  // rho gt 1
  inp << 0.3, 0.3, 1.3;
  invalid_args.push_back(inp);
  test_invalid_args(invalid_args);
}
TEST(MathFunctions, binormal_copula_cdf_throw_bad_sizes_all_types) {
  Matrix<double, Dynamic, 1> inp(4);
  vector<Matrix<double, Dynamic, 1>> invalid_args;
  // RV1 nan
  inp << 0.1, 0.2, 0.4, 0.5;
  invalid_args.push_back(inp);
  inp.resize(2);
  inp << 0.1, 0.5;
  invalid_args.push_back(inp);
  test_bad_sizes(invalid_args);
}
TEST(MathFunctions, binormal_copula_cdf_no_throw_valid_args_all_types) {
  Matrix<double, Dynamic, 1> inp(3);
  vector<Matrix<double, Dynamic, 1>> valid_args;
  inp << 0.1, 0.7, 0.3;
  valid_args.push_back(inp);
  inp << 0.3, 0.1, 0.3;
  valid_args.push_back(inp);
  inp << 0.3, 0.3, 0.1;
  valid_args.push_back(inp);
  inp << 0.9, 0.3, 0.3;
  valid_args.push_back(inp);
  inp << 0.3, 0.3, 0.9;
  valid_args.push_back(inp);
  inp << 1, 0.3, 0.9;
  valid_args.push_back(inp);
  inp << 0, 0.3, 0.9;
  valid_args.push_back(inp);
  inp << 0.3, 1, 0.9;
  valid_args.push_back(inp);
  inp << 0.3, 0, 0.9;
  valid_args.push_back(inp);
  test_valid_args(valid_args);
}

TEST(MathFunctions, binormal_copula_cdf_boundaries_test) {
  using std::asin;
  using std::exp;
  using std::log;
  Matrix<double, Dynamic, 1> inp(3);
  vector<Matrix<double, Dynamic, 1>> valid_args;
  vector<double> vals;

  // Independent RVs
  inp << 0.1, 0.7, 0.0;
  valid_args.push_back(inp);
  vals.push_back(0.1 * 0.7);

  // Perfectly correlated RVs
  inp << 0.1, 0.7, 1.0;
  valid_args.push_back(inp);
  vals.push_back(0.1);

  // Perfectly anticorrelated RVs
  inp << 0.5, 0.7, -1.0;
  valid_args.push_back(inp);
  vals.push_back(0.5 + 0.7 - 1.0);

  // Perfectly anticorrelated RVs
  inp << 0.1, 0.7, -1.0;
  valid_args.push_back(inp);
  vals.push_back(0);

  // Integrate out first RV
  inp << 1, 0.7, 0.5;
  valid_args.push_back(inp);
  vals.push_back(0.7);

  // Integrate out second RV
  inp << 0.7, 1, 0.5;
  valid_args.push_back(inp);
  vals.push_back(0.7);

  // Zero out CDF
  inp << 0, 0.7, 0.5;
  valid_args.push_back(inp);
  vals.push_back(0);

  // Zero out CDF
  inp << 0.7, 0, 0.5;
  valid_args.push_back(inp);
  vals.push_back(0);

  test_vals(valid_args, vals);
}
TEST(MathFunctions, vec_binormal_integral_val_and_grad_test_var_var) {
  V_R_binorm_copula_cdf<Dynamic, 1> dist_fun_V_R;
  V_R_binorm_copula_cdf<1, Dynamic> dist_fun_RV_R;
  V_R_binorm_copula_manual tru_fun;
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << 0.4, 0.7, 0.8;
  compare_grad(dist_fun_V_R, tru_fun, inp_vec);
  compare_fd(dist_fun_V_R, inp_vec);

  compare_grad(dist_fun_RV_R, tru_fun, inp_vec);
  compare_fd(dist_fun_RV_R, inp_vec);
}
TEST(MathFunctions, vec_binormal_integral_val_and_grad_test_var_doub) {
  V_D_binorm_copula_cdf<Dynamic, 1> dist_fun_V_R(0.3);
  V_D_binorm_copula_cdf<1, Dynamic> dist_fun_RV_R(0.3);
  V_D_binorm_copula_manual tru_fun(0.3);
  Matrix<double, Dynamic, 1> inp_vec(2);
  inp_vec << 0.4, 0.7;
  compare_grad(dist_fun_V_R, tru_fun, inp_vec);
  compare_fd(dist_fun_V_R, inp_vec);

  compare_grad(dist_fun_RV_R, tru_fun, inp_vec);
  compare_fd(dist_fun_RV_R, inp_vec);
}
TEST(MathFunctions, vec_binormal_integral_val_and_grad_test_var_doub_bound) {
  using stan::math::inv_Phi;
  using std::exp;
  using std::pow;
  using std::sqrt;
  V_D_binorm_copula_cdf<Dynamic, 1> dist_fun_V_R(0.3);
  V_D_binorm_copula_cdf<Dynamic, 1> dist_fun_V_R_rho_1(1.0);
  V_D_binorm_copula_cdf<Dynamic, 1> dist_fun_V_R_rho_neg_1(-1.0);
  V_D_binorm_copula_cdf<1, Dynamic> dist_fun_RV_R(0.3);
  Matrix<double, Dynamic, 1> inp_vec(2);
  inp_vec << 0, 0;
  Matrix<double, Dynamic, 1> expected_grad(2);
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R, expected_grad, inp_vec);

  inp_vec << 0, 0.3;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R, expected_grad, inp_vec);

  inp_vec << 0.3, 0;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R, expected_grad, inp_vec);

  inp_vec << 1, 1;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R, expected_grad, inp_vec);

  inp_vec << 1, 0.3;
  expected_grad << 0, 1;
  compare_grad_known(dist_fun_V_R, expected_grad, inp_vec);

  inp_vec << 0.3, 1;
  expected_grad << 1, 0;
  compare_grad_known(dist_fun_V_R, expected_grad, inp_vec);

  inp_vec << 1, 0;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R, expected_grad, inp_vec);

  inp_vec << 0, 1;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R, expected_grad, inp_vec);

  // rho = 1
  inp_vec << 0, 0;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);

  inp_vec << 0, 1;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);

  inp_vec << 1, 0;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);

  inp_vec << 1, 1;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);

  inp_vec << 1, 0.3;
  expected_grad << 0, 1;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);

  inp_vec << 0.3, 1;
  expected_grad << 1, 0;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);

  inp_vec << 0.3, 0;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);

  inp_vec << 0, 0.3;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);

  inp_vec << 0.4, 0.3;
  expected_grad << 0, 1;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);

  inp_vec << 0.3, 0.4;
  expected_grad << 1, 0;
  compare_grad_known(dist_fun_V_R_rho_1, expected_grad, inp_vec);
  //
  // rho = -1
  inp_vec << 0, 0;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);

  inp_vec << 0, 1;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);

  inp_vec << 1, 0;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);

  inp_vec << 1, 1;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);

  inp_vec << 1, 0.3;
  expected_grad << 0, 1;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);

  inp_vec << 0.3, 1;
  expected_grad << 1, 0;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);

  inp_vec << 0.3, 0;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);

  inp_vec << 0, 0.3;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);

  inp_vec << 0.4, 0.7;
  expected_grad << 1, 1;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);

  inp_vec << 0.3, 0.4;
  expected_grad << 0, 0;
  compare_grad_known(dist_fun_V_R_rho_neg_1, expected_grad, inp_vec);
}
TEST(MathFunctions, vec_binormal_integral_val_and_grad_test_doub_var) {
  Matrix<double, Dynamic, 1> y1(2);
  y1(0) = 0.2;
  y1(1) = 0.9;
  Matrix<double, 1, Dynamic> y2 = y1.transpose();
  VD_R_binorm_copula_cdf<Dynamic, 1> dist_fun_V_R(y1);
  VD_R_binorm_copula_cdf<1, Dynamic> dist_fun_RV_R(y2);
  VD_R_binorm_copula_manual<Dynamic, 1> tru_fun(y1);
  Matrix<double, Dynamic, 1> inp_vec(1);
  inp_vec(0) = 0.9;
  compare_grad(dist_fun_V_R, tru_fun, inp_vec);
  compare_fd(dist_fun_V_R, inp_vec);

  compare_grad(dist_fun_RV_R, tru_fun, inp_vec);
  compare_fd(dist_fun_RV_R, inp_vec);
}
