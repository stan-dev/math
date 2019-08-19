#include <stan/math/prim/scal.hpp>
#include <Eigen/Core>
#include <gtest/gtest.h>
#include <vector>
#include <type_traits>

// for testing that stan with just prim/scal does not know what fvar is.
template <typename T>
class fvar {};

// dummy var class
class var {};

// Note: Failure to find specialization is defined as "false"
template <typename T, typename = void>
struct disable_if_tester : std::false_type {};

template <typename T>
struct disable_if_tester<T, stan::disable_if<T>> : std::true_type {};

TEST(template_enablers, disable_if) {
  auto val = disable_if_tester<std::false_type>::value;
  EXPECT_TRUE(val);
  val = disable_if_tester<std::true_type>::value;
  EXPECT_FALSE(val);
}

template <typename T1, typename T2, typename = void>
struct disable_if_all_tester : std::false_type {};

template <typename T1, typename T2>
struct disable_if_all_tester<T1, T2, stan::disable_if_all<T1, T2>>
    : std::true_type {};

TEST(template_enablers, disable_if_all) {
  auto val = disable_if_all_tester<std::false_type, std::false_type>::value;
  EXPECT_TRUE(val);
  val = disable_if_all_tester<std::true_type, std::true_type>::value;
  EXPECT_FALSE(val);
  val = disable_if_all_tester<std::true_type, std::false_type>::value;
  EXPECT_TRUE(val);
}

template <typename T1, typename T2, typename = void>
struct disable_if_any_tester : std::false_type {};

template <typename T1, typename T2>
struct disable_if_any_tester<T1, T2, stan::disable_if_any<T1, T2>>
    : std::true_type {};

TEST(template_enablers, disable_if_any) {
  auto val = disable_if_any_tester<std::false_type, std::false_type>::value;
  EXPECT_TRUE(val);
  val = disable_if_any_tester<std::true_type, std::true_type>::value;
  EXPECT_FALSE(val);
  val = disable_if_any_tester<std::true_type, std::false_type>::value;
  EXPECT_FALSE(val);
}

template <typename T1, typename T2, typename = void>
struct enable_if_all_tester : std::false_type {};

template <typename T1, typename T2>
struct enable_if_all_tester<T1, T2, stan::enable_if_all<T1, T2>>
    : std::true_type {};

TEST(template_enablers, enable_if_all) {
  auto val = enable_if_all_tester<std::false_type, std::false_type>::value;
  EXPECT_FALSE(val);
  val = enable_if_all_tester<std::true_type, std::true_type>::value;
  EXPECT_TRUE(val);
  val = enable_if_all_tester<std::true_type, std::false_type>::value;
  EXPECT_FALSE(val);
}

template <typename T1, typename T2, typename = void>
struct enable_if_any_tester : std::false_type {};

template <typename T1, typename T2>
struct enable_if_any_tester<T1, T2, stan::enable_if_any<T1, T2>>
    : std::true_type {};

TEST(template_enablers, enable_if_any) {
  auto val = enable_if_any_tester<std::false_type, std::false_type>::value;
  EXPECT_FALSE(val);
  val = enable_if_any_tester<std::true_type, std::true_type>::value;
  EXPECT_TRUE(val);
  val = enable_if_any_tester<std::true_type, std::false_type>::value;
  EXPECT_TRUE(val);
}

TEST(is_check, is_contains_arithmetic) {
  using stan::is_arithmetic_container;
  class var {};
  auto val = is_arithmetic_container<double>::value;
  EXPECT_TRUE(val);

  val = is_arithmetic_container<var>::value;
  EXPECT_FALSE(val);

  val = is_arithmetic_container<std::vector<double>>::value;
  EXPECT_TRUE(val);

  val = is_arithmetic_container<Eigen::Matrix<double, -1, 1>>::value;
  EXPECT_TRUE(val);
}

TEST(is_check, is_contains_floating_point) {
  using stan::is_floating_point_container;
  auto val = is_floating_point_container<double>::value;
  EXPECT_TRUE(val);

  val = is_floating_point_container<int>::value;
  EXPECT_FALSE(val);

  val = is_floating_point_container<std::vector<double>>::value;
  EXPECT_TRUE(val);

  val = is_floating_point_container<Eigen::Matrix<double, -1, 1>>::value;
  EXPECT_TRUE(val);
}

TEST(is_check, is_contains_var) {
  using stan::is_var_container;

  auto val = is_var_container<double>::value;
  EXPECT_FALSE(val);
  val = is_var_container<var>::value;
  EXPECT_FALSE(val);
  val = is_var_container<fvar<double>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<fvar<var>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<fvar<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<fvar<fvar<var>>>::value;
  EXPECT_FALSE(val);

  val = is_var_container<std::vector<double>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<std::vector<var>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<std::vector<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<std::vector<fvar<var>>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<std::vector<fvar<fvar<double>>>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<std::vector<fvar<fvar<var>>>>::value;
  EXPECT_FALSE(val);

  val = is_var_container<Eigen::Matrix<double, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<Eigen::Matrix<var, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<Eigen::Matrix<fvar<double>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<Eigen::Matrix<fvar<var>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<Eigen::Matrix<fvar<fvar<double>>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_var_container<Eigen::Matrix<fvar<fvar<var>>, -1, -1>>::value;
  EXPECT_FALSE(val);
}

TEST(is_check, is_contains_fvar) {
  using stan::is_fvar_container;

  auto val = is_fvar_container<double>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<var>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<fvar<double>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<fvar<var>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<fvar<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<fvar<fvar<var>>>::value;
  EXPECT_FALSE(val);

  val = is_fvar_container<std::vector<double>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<std::vector<var>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<std::vector<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<std::vector<fvar<var>>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<std::vector<fvar<fvar<double>>>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<std::vector<fvar<fvar<var>>>>::value;
  EXPECT_FALSE(val);

  val = is_fvar_container<Eigen::Matrix<double, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<Eigen::Matrix<var, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<Eigen::Matrix<fvar<double>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<Eigen::Matrix<fvar<var>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<Eigen::Matrix<fvar<fvar<double>>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_fvar_container<Eigen::Matrix<fvar<fvar<var>>, -1, -1>>::value;
  EXPECT_FALSE(val);
}

TEST(is_check, is_contains_ad_type) {
  using stan::is_ad_type_container;

  auto val = is_ad_type_container<double>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<var>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<fvar<double>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<fvar<var>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<fvar<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<fvar<fvar<var>>>::value;
  EXPECT_FALSE(val);

  val = is_ad_type_container<std::vector<double>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<std::vector<var>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<std::vector<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<std::vector<fvar<var>>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<std::vector<fvar<fvar<double>>>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<std::vector<fvar<fvar<var>>>>::value;
  EXPECT_FALSE(val);

  val = is_ad_type_container<Eigen::Matrix<double, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<Eigen::Matrix<var, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<Eigen::Matrix<fvar<double>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<Eigen::Matrix<fvar<var>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<Eigen::Matrix<fvar<fvar<double>>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<Eigen::Matrix<fvar<fvar<var>>, -1, -1>>::value;
  EXPECT_FALSE(val);
}

TEST(is_check, is_eigen_arithmetic) {
  using stan::is_eigen_arithmetic;

  auto val = is_eigen_arithmetic<double>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<var>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<fvar<double>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<fvar<var>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<fvar<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<fvar<fvar<var>>>::value;
  EXPECT_FALSE(val);

  val = is_eigen_arithmetic<std::vector<double>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<std::vector<var>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<std::vector<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<std::vector<fvar<var>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<std::vector<fvar<fvar<double>>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<std::vector<fvar<fvar<var>>>>::value;
  EXPECT_FALSE(val);

  val = is_eigen_arithmetic<Eigen::Matrix<double, -1, -1>>::value;
  EXPECT_TRUE(val);
  val = is_eigen_arithmetic<Eigen::Matrix<var, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<Eigen::Matrix<fvar<double>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<Eigen::Matrix<fvar<var>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<Eigen::Matrix<fvar<fvar<double>>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_arithmetic<Eigen::Matrix<fvar<fvar<var>>, -1, -1>>::value;
  EXPECT_FALSE(val);
}

TEST(is_check, is_eigen_var) {
  using stan::is_eigen_var;

  auto val = is_eigen_var<double>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<var>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<fvar<double>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<fvar<var>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<fvar<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<fvar<fvar<var>>>::value;
  EXPECT_FALSE(val);

  val = is_eigen_var<std::vector<double>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<std::vector<var>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<std::vector<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<std::vector<fvar<var>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<std::vector<fvar<fvar<double>>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<std::vector<fvar<fvar<var>>>>::value;
  EXPECT_FALSE(val);

  val = is_eigen_var<Eigen::Matrix<double, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<Eigen::Matrix<var, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<Eigen::Matrix<fvar<double>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<Eigen::Matrix<fvar<var>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<Eigen::Matrix<fvar<fvar<double>>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_var<Eigen::Matrix<fvar<fvar<var>>, -1, -1>>::value;
  EXPECT_FALSE(val);
}

TEST(is_check, is_eigen_fvar) {
  using stan::is_eigen_fvar;

  auto val = is_eigen_fvar<double>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<var>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<fvar<double>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<fvar<var>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<fvar<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<fvar<fvar<var>>>::value;
  EXPECT_FALSE(val);

  val = is_eigen_fvar<std::vector<double>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<std::vector<var>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<std::vector<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<std::vector<fvar<var>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<std::vector<fvar<fvar<double>>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<std::vector<fvar<fvar<var>>>>::value;
  EXPECT_FALSE(val);

  val = is_eigen_fvar<Eigen::Matrix<double, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<Eigen::Matrix<var, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<Eigen::Matrix<fvar<double>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<Eigen::Matrix<fvar<var>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<Eigen::Matrix<fvar<fvar<double>>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_fvar<Eigen::Matrix<fvar<fvar<var>>, -1, -1>>::value;
  EXPECT_FALSE(val);
}

TEST(is_check, is_contains_var_or_arithmetic) {
  using stan::is_var_or_arithmetic_container;
  class var {};
  auto val = is_var_or_arithmetic_container<var>::value;
  EXPECT_FALSE(val);

  val = is_var_or_arithmetic_container<int>::value;
  EXPECT_TRUE(val);

  val = is_var_or_arithmetic_container<std::vector<double>>::value;
  EXPECT_TRUE(val);

  val = is_var_or_arithmetic_container<std::vector<double>>::value;
  EXPECT_TRUE(val);

  val = is_var_or_arithmetic_container<Eigen::Matrix<double, -1, 1>>::value;
  EXPECT_TRUE(val);
}

TEST(is_check, is_contains_stan_scalar) {
  using stan::is_stan_scalar_container;
  class var {};

  auto val = is_stan_scalar_container<var>::value;
  EXPECT_FALSE(val);

  val = is_stan_scalar_container<double>::value;
  EXPECT_TRUE(val);

  val = is_stan_scalar_container<std::vector<double>>::value;
  EXPECT_TRUE(val);

  val = is_stan_scalar_container<std::vector<fvar<double>>>::value;
  EXPECT_FALSE(val);

  val = is_stan_scalar_container<Eigen::Matrix<double, -1, 1>>::value;
  EXPECT_TRUE(val);
}

TEST(is_check, is_eigen_or_stan_scalar) {
  using stan::is_eigen_or_stan_scalar;
  auto val = is_eigen_or_stan_scalar<double>::value;
  EXPECT_TRUE(val);
  val = is_eigen_or_stan_scalar<var>::value;
  EXPECT_FALSE(val);
  val = is_eigen_or_stan_scalar<fvar<double>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_or_stan_scalar<fvar<var>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_or_stan_scalar<fvar<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_or_stan_scalar<fvar<fvar<var>>>::value;
  EXPECT_FALSE(val);

  val = is_eigen_or_stan_scalar<std::vector<double>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_or_stan_scalar<std::vector<var>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_or_stan_scalar<std::vector<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_or_stan_scalar<std::vector<fvar<var>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_or_stan_scalar<std::vector<fvar<fvar<double>>>>::value;
  EXPECT_FALSE(val);
  val = is_eigen_or_stan_scalar<std::vector<fvar<fvar<var>>>>::value;
  EXPECT_FALSE(val);

  val = is_eigen_or_stan_scalar<Eigen::Matrix<double, -1, -1>>::value;
  EXPECT_TRUE(val);
  val = is_eigen_or_stan_scalar<Eigen::Matrix<var, -1, -1>>::value;
  EXPECT_TRUE(val);
  val = is_eigen_or_stan_scalar<Eigen::Matrix<fvar<double>, -1, -1>>::value;
  EXPECT_TRUE(val);
  val = is_eigen_or_stan_scalar<Eigen::Matrix<fvar<var>, -1, -1>>::value;
  EXPECT_TRUE(val);
  val = is_eigen_or_stan_scalar<
      Eigen::Matrix<fvar<fvar<double>>, -1, -1>>::value;
  EXPECT_TRUE(val);
  val = is_eigen_or_stan_scalar<Eigen::Matrix<fvar<fvar<var>>, -1, -1>>::value;
  EXPECT_TRUE(val);
}

TEST(is_check, is_dot_product) {
  using stan::is_dot_product;

  using left_dot = typename Eigen::Matrix<double, 1, -1>;
  using right_dot = typename Eigen::Matrix<double, -1, 1>;
  auto val = is_dot_product<left_dot, right_dot>::value;
  EXPECT_TRUE(val);

  val = is_dot_product<right_dot, left_dot>::value;
  EXPECT_FALSE(val);

  val = is_dot_product<double, right_dot>::value;
  EXPECT_FALSE(val);

  using full_mat = typename Eigen::Matrix<double, -1, -1>;
  val = is_dot_product<full_mat, right_dot>::value;
  EXPECT_FALSE(val);
}

TEST(is_check, is_eigen_rows_match) {
  using stan::is_eigen_rows_match;

  using left_dot = typename Eigen::Matrix<double, 1, -1>;
  using right_dot = typename Eigen::Matrix<double, -1, 1>;
  auto val = is_eigen_rows_match<left_dot, right_dot>::value;
  EXPECT_FALSE(val);

  val = is_eigen_rows_match<double, right_dot>::value;
  EXPECT_FALSE(val);

  using full_mat = typename Eigen::Matrix<double, -1, -1>;
  val = is_eigen_rows_match<full_mat, right_dot>::value;
  EXPECT_TRUE(val);
}

TEST(is_check, is_eigen_cols_match) {
  using stan::is_eigen_cols_match;

  using left_dot = typename Eigen::Matrix<double, 1, -1>;
  using right_dot = typename Eigen::Matrix<double, -1, 1>;
  auto val = is_eigen_cols_match<left_dot, right_dot>::value;
  EXPECT_FALSE(val);

  val = is_eigen_cols_match<double, right_dot>::value;
  EXPECT_FALSE(val);

  using full_mat = typename Eigen::Matrix<double, -1, -1>;
  val = is_eigen_cols_match<full_mat, right_dot>::value;
  EXPECT_FALSE(val);

  val = is_eigen_cols_match<full_mat, left_dot>::value;
  EXPECT_TRUE(val);
}

TEST(is_check, is_eigen_dims_match) {
  using stan::is_eigen_dims_match;

  using left_dot = typename Eigen::Matrix<double, 1, -1>;
  using right_dot = typename Eigen::Matrix<double, -1, 1>;
  auto val = is_eigen_dims_match<left_dot, right_dot>::value;
  EXPECT_FALSE(val);

  val = is_eigen_dims_match<double, right_dot>::value;
  EXPECT_FALSE(val);

  using full_mat = typename Eigen::Matrix<double, -1, -1>;
  val = is_eigen_dims_match<full_mat, full_mat>::value;
  EXPECT_TRUE(val);

  val = is_eigen_dims_match<left_dot, left_dot>::value;
  EXPECT_TRUE(val);
}
