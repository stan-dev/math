#include <test/unit/math/prim/mat/prob/binormal_copula_cdf_test_helper.hpp>
#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::fvar;
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
    EXPECT_FLOAT_EQ(grad_f1(i), grad_f2(i));
}

template <class F1, class F2>
void compare_hess(const F1& f1, const F2& f2,
                  const Matrix<double, Dynamic, 1>& inp_vec) {
  double f1_eval;
  Matrix<double, Dynamic, 1> grad_f1;
  Matrix<double, Dynamic, Dynamic> H_f1;
  double f2_eval;
  Matrix<double, Dynamic, 1> grad_f2;
  Matrix<double, Dynamic, Dynamic> H_f2;
  stan::math::hessian(f1, inp_vec, f1_eval, grad_f1, H_f1);
  stan::math::hessian(f2, inp_vec, f2_eval, grad_f2, H_f2);
  EXPECT_FLOAT_EQ(f1_eval, f2_eval);
  for (int i = 0; i < grad_f1.size(); ++i)
    EXPECT_FLOAT_EQ(grad_f1(i), grad_f2(i));
  for (int i = 0; i < H_f1.rows(); ++i)
    for (int j = 0; j < H_f1.cols(); ++j)
      EXPECT_FLOAT_EQ(H_f1(i, j), H_f2(i, j));
}

void test_invalid_args(const vector<Matrix<double, Dynamic, 1>>& invalid_args) {
  for (size_t i = 0; i < invalid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = invalid_args[i];
    auto args_V_ff = make_args<fvar<double>, fvar<double>, Dynamic, 1>(inp);
    auto args_V_df = make_args<double, fvar<double>, Dynamic, 1>(inp);
    auto args_V_fd = make_args<fvar<double>, double, Dynamic, 1>(inp);
    auto args_RV_ff = make_args<fvar<double>, fvar<double>, 1, Dynamic>(inp);
    auto args_RV_df = make_args<double, fvar<double>, 1, Dynamic>(inp);
    auto args_RV_fd = make_args<fvar<double>, double, 1, Dynamic>(inp);
    auto args_V_ffff
        = make_args<fvar<fvar<double>>, fvar<fvar<double>>, Dynamic, 1>(inp);
    auto args_V_dff = make_args<double, fvar<fvar<double>>, Dynamic, 1>(inp);
    auto args_V_ffd = make_args<fvar<fvar<double>>, double, Dynamic, 1>(inp);
    auto args_RV_ffff
        = make_args<fvar<fvar<double>>, fvar<fvar<double>>, 1, Dynamic>(inp);
    auto args_RV_dff = make_args<double, fvar<fvar<double>>, 1, Dynamic>(inp);
    auto args_RV_ffd = make_args<fvar<fvar<double>>, double, 1, Dynamic>(inp);
    EXPECT_THROW(call_args(args_V_ff), std::domain_error);
    EXPECT_THROW(call_args(args_V_df), std::domain_error);
    EXPECT_THROW(call_args(args_V_fd), std::domain_error);
    EXPECT_THROW(call_args(args_RV_ff), std::domain_error);
    EXPECT_THROW(call_args(args_RV_df), std::domain_error);
    EXPECT_THROW(call_args(args_RV_fd), std::domain_error);
    EXPECT_THROW(call_args(args_V_ffff), std::domain_error);
    EXPECT_THROW(call_args(args_V_dff), std::domain_error);
    EXPECT_THROW(call_args(args_V_ffd), std::domain_error);
    EXPECT_THROW(call_args(args_RV_ffff), std::domain_error);
    EXPECT_THROW(call_args(args_RV_dff), std::domain_error);
    EXPECT_THROW(call_args(args_RV_ffd), std::domain_error);
  }
}

void test_bad_sizes(const vector<Matrix<double, Dynamic, 1>>& invalid_args) {
  for (size_t i = 0; i < invalid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = invalid_args[i];
    auto args_V_ff = make_args<fvar<double>, fvar<double>, Dynamic, 1>(inp);
    auto args_V_df = make_args<double, fvar<double>, Dynamic, 1>(inp);
    auto args_V_fd = make_args<fvar<double>, double, Dynamic, 1>(inp);
    auto args_RV_ff = make_args<fvar<double>, fvar<double>, 1, Dynamic>(inp);
    auto args_RV_df = make_args<double, fvar<double>, 1, Dynamic>(inp);
    auto args_RV_fd = make_args<fvar<double>, double, 1, Dynamic>(inp);
    auto args_V_ffff
        = make_args<fvar<fvar<double>>, fvar<fvar<double>>, Dynamic, 1>(inp);
    auto args_V_dff = make_args<double, fvar<fvar<double>>, Dynamic, 1>(inp);
    auto args_V_ffd = make_args<fvar<fvar<double>>, double, Dynamic, 1>(inp);
    auto args_RV_ffff
        = make_args<fvar<fvar<double>>, fvar<fvar<double>>, 1, Dynamic>(inp);
    auto args_RV_dff = make_args<double, fvar<fvar<double>>, 1, Dynamic>(inp);
    auto args_RV_ffd = make_args<fvar<fvar<double>>, double, 1, Dynamic>(inp);
    EXPECT_THROW(call_args(args_V_ff), std::invalid_argument);
    EXPECT_THROW(call_args(args_V_df), std::invalid_argument);
    EXPECT_THROW(call_args(args_V_fd), std::invalid_argument);
    EXPECT_THROW(call_args(args_RV_ff), std::invalid_argument);
    EXPECT_THROW(call_args(args_RV_df), std::invalid_argument);
    EXPECT_THROW(call_args(args_RV_fd), std::invalid_argument);
    EXPECT_THROW(call_args(args_V_ffff), std::invalid_argument);
    EXPECT_THROW(call_args(args_V_dff), std::invalid_argument);
    EXPECT_THROW(call_args(args_V_ffd), std::invalid_argument);
    EXPECT_THROW(call_args(args_RV_ffff), std::invalid_argument);
    EXPECT_THROW(call_args(args_RV_dff), std::invalid_argument);
    EXPECT_THROW(call_args(args_RV_ffd), std::invalid_argument);
  }
}

void test_valid_args(const vector<Matrix<double, Dynamic, 1>>& valid_args) {
  for (size_t i = 0; i < valid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = valid_args[i];
    auto args_V_ff = make_args<fvar<double>, fvar<double>, Dynamic, 1>(inp);
    auto args_V_df = make_args<double, fvar<double>, Dynamic, 1>(inp);
    auto args_V_fd = make_args<fvar<double>, double, Dynamic, 1>(inp);
    auto args_RV_ff = make_args<fvar<double>, fvar<double>, 1, Dynamic>(inp);
    auto args_RV_df = make_args<double, fvar<double>, 1, Dynamic>(inp);
    auto args_RV_fd = make_args<fvar<double>, double, 1, Dynamic>(inp);
    auto args_V_ffff
        = make_args<fvar<fvar<double>>, fvar<fvar<double>>, Dynamic, 1>(inp);
    auto args_V_dff = make_args<double, fvar<fvar<double>>, Dynamic, 1>(inp);
    auto args_V_ffd = make_args<fvar<fvar<double>>, double, Dynamic, 1>(inp);
    auto args_RV_ffff
        = make_args<fvar<fvar<double>>, fvar<fvar<double>>, 1, Dynamic>(inp);
    auto args_RV_dff = make_args<double, fvar<fvar<double>>, 1, Dynamic>(inp);
    auto args_RV_ffd = make_args<fvar<fvar<double>>, double, 1, Dynamic>(inp);
    EXPECT_NO_THROW(call_args(args_V_ff));
    EXPECT_NO_THROW(call_args(args_V_df));
    EXPECT_NO_THROW(call_args(args_V_fd));
    EXPECT_NO_THROW(call_args(args_RV_ff));
    EXPECT_NO_THROW(call_args(args_RV_df));
    EXPECT_NO_THROW(call_args(args_RV_fd));
    EXPECT_NO_THROW(call_args(args_V_ffff));
    EXPECT_NO_THROW(call_args(args_V_dff));
    EXPECT_NO_THROW(call_args(args_V_ffd));
    EXPECT_NO_THROW(call_args(args_RV_ffff));
    EXPECT_NO_THROW(call_args(args_RV_dff));
    EXPECT_NO_THROW(call_args(args_RV_ffd));
  }
}

void test_vals(const vector<Matrix<double, Dynamic, 1>>& valid_args,
               const vector<double>& vals) {
  for (size_t i = 0; i < valid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = valid_args[i];
    auto args_V_ff = make_args<fvar<double>, fvar<double>, Dynamic, 1>(inp);
    auto args_V_df = make_args<double, fvar<double>, Dynamic, 1>(inp);
    auto args_V_fd = make_args<fvar<double>, double, Dynamic, 1>(inp);
    auto args_RV_ff = make_args<fvar<double>, fvar<double>, 1, Dynamic>(inp);
    auto args_RV_df = make_args<double, fvar<double>, 1, Dynamic>(inp);
    auto args_RV_fd = make_args<fvar<double>, double, 1, Dynamic>(inp);
    auto args_V_ffff
        = make_args<fvar<fvar<double>>, fvar<fvar<double>>, Dynamic, 1>(inp);
    auto args_V_dff = make_args<double, fvar<fvar<double>>, Dynamic, 1>(inp);
    auto args_V_ffd = make_args<fvar<fvar<double>>, double, Dynamic, 1>(inp);
    auto args_RV_ffff
        = make_args<fvar<fvar<double>>, fvar<fvar<double>>, 1, Dynamic>(inp);
    auto args_RV_dff = make_args<double, fvar<fvar<double>>, 1, Dynamic>(inp);
    auto args_RV_ffd = make_args<fvar<fvar<double>>, double, 1, Dynamic>(inp);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_V_ff)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_V_df)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_V_fd)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_RV_ff)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_RV_df)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_RV_fd)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_V_ffff)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_V_dff)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_V_ffd)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_RV_ffff)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_RV_dff)), vals[i]);
    EXPECT_FLOAT_EQ(stan::math::value_of_rec(call_args(args_RV_ffd)), vals[i]);
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
TEST(MathFunctions, vec_binormal_integral_val_and_grad_test_fvar_fvar) {
  V_R_binorm_copula_cdf<Dynamic, 1> dist_fun_V_R;
  V_R_binorm_copula_cdf<1, Dynamic> dist_fun_RV_R;
  V_R_binorm_copula_manual tru_fun;
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << 0.4, 0.7, 0.8;
  compare_grad(dist_fun_V_R, tru_fun, inp_vec);

  compare_grad(dist_fun_RV_R, tru_fun, inp_vec);

  compare_hess(dist_fun_RV_R, tru_fun, inp_vec);
  compare_hess(dist_fun_V_R, tru_fun, inp_vec);
}
TEST(MathFunctions, vec_binormal_integral_val_and_grad_test_fvar_doub) {
  V_D_binorm_copula_cdf<Dynamic, 1> dist_fun_V_R(0.3);
  V_D_binorm_copula_cdf<1, Dynamic> dist_fun_RV_R(0.3);
  V_D_binorm_copula_manual tru_fun(0.3);
  Matrix<double, Dynamic, 1> inp_vec(2);
  inp_vec << 0.4, 0.7;
  compare_grad(dist_fun_V_R, tru_fun, inp_vec);

  compare_grad(dist_fun_RV_R, tru_fun, inp_vec);

  compare_hess(dist_fun_RV_R, tru_fun, inp_vec);
  compare_hess(dist_fun_V_R, tru_fun, inp_vec);
}
TEST(MathFunctions, vec_binormal_integral_val_and_grad_test_doub_fvar) {
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

  compare_grad(dist_fun_RV_R, tru_fun, inp_vec);

  compare_hess(dist_fun_RV_R, tru_fun, inp_vec);
  compare_hess(dist_fun_V_R, tru_fun, inp_vec);
}
