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
void compare_vals(const F1& f1, const F2& f2,
                  const Matrix<double, Dynamic, 1>& inp_vec) {
  double f1_eval = f1(inp_vec);
  double f2_eval = f2(inp_vec);
  EXPECT_FLOAT_EQ(f1_eval, f2_eval);
}

void test_invalid_args(const vector<Matrix<double, Dynamic, 1>>& invalid_args) {
  for (size_t i = 0; i < invalid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = invalid_args[i];
    auto args_V_dd = make_args<double, double, Dynamic, 1>(inp);
    auto args_RV_dd = make_args<double, double, 1, Dynamic>(inp);
    EXPECT_THROW(call_args(args_V_dd), std::domain_error);
    EXPECT_THROW(call_args(args_RV_dd), std::domain_error);
  }
}

void test_valid_args(const vector<Matrix<double, Dynamic, 1>>& valid_args) {
  for (size_t i = 0; i < valid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = valid_args[i];
    auto args_V_dd = make_args<double, double, Dynamic, 1>(inp);
    auto args_RV_dd = make_args<double, double, 1, Dynamic>(inp);
    EXPECT_NO_THROW(call_args(args_V_dd));
    EXPECT_NO_THROW(call_args(args_RV_dd));
  }
}

void test_vals(const vector<Matrix<double, Dynamic, 1>>& valid_args,
               const vector<double>& vals) {
  for (size_t i = 0; i < valid_args.size(); ++i) {
    Matrix<double, Dynamic, 1> inp = valid_args[i];
    auto args_V_dd = make_args<double, double, Dynamic, 1>(inp);
    auto args_RV_dd = make_args<double, double, 1, Dynamic>(inp);
    EXPECT_FLOAT_EQ(call_args(args_V_dd), vals[i]);
    EXPECT_FLOAT_EQ(call_args(args_RV_dd), vals[i]);
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
TEST(MathFunctions, binormal_copula_cdf_size_test) {
  Matrix<double, Dynamic, 1> inv_y_V(3);
  Matrix<double, 1, Dynamic> inv_y_RV(3);
  double rho = 0.3;
  inv_y_V.setZero();
  inv_y_RV.setZero();
  EXPECT_THROW(stan::math::binormal_copula_cdf(inv_y_V, rho),
               std::invalid_argument);
  EXPECT_THROW(stan::math::binormal_copula_cdf(inv_y_RV, rho),
               std::invalid_argument);
  inv_y_V.resize(1);
  inv_y_RV.resize(1);
  inv_y_V.setZero();
  inv_y_RV.setZero();
  EXPECT_THROW(stan::math::binormal_copula_cdf(inv_y_V, rho),
               std::invalid_argument);
  EXPECT_THROW(stan::math::binormal_copula_cdf(inv_y_RV, rho),
               std::invalid_argument);
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
TEST(MathFunctions, vec_binormal_integral_val_test) {
  V_R_binorm_copula_cdf<Dynamic, 1> dist_fun_V_R;
  V_R_binorm_copula_cdf<1, Dynamic> dist_fun_RV_R;
  V_R_binorm_copula_manual tru_fun;
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << 0.4, 0.7, 0.8;
  compare_vals(dist_fun_V_R, tru_fun, inp_vec);
  compare_vals(dist_fun_RV_R, tru_fun, inp_vec);
}
