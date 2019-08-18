#include <stan/math/rev/scal.hpp>
#include <Eigen/Core>
#include <gtest/gtest.h>
#include <vector>
#include <type_traits>

// These are the is_* methods that check for var types

// for testing that stan with just prim/scal does not know what fvar is.
template <typename T>
class fvar {};

TEST(is_check, is_contains_var) {
  using stan::is_var_container;
  using stan::math::var;

  auto val = is_var_container<double>::value;
  EXPECT_FALSE(val);
  val = is_var_container<var>::value;
  EXPECT_TRUE(val);
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
  EXPECT_TRUE(val);
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
  EXPECT_TRUE(val);
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
  using stan::math::var;

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

TEST(is_check, is_contains_stan_scalar) {
  using stan::is_stan_scalar_container;
  using stan::math::var;

  auto val = is_stan_scalar_container<double>::value;
  EXPECT_TRUE(val);
  val = is_stan_scalar_container<var>::value;
  EXPECT_TRUE(val);
  val = is_stan_scalar_container<fvar<double>>::value;
  EXPECT_FALSE(val);
  val = is_stan_scalar_container<fvar<var>>::value;
  EXPECT_FALSE(val);
  val = is_stan_scalar_container<fvar<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_stan_scalar_container<fvar<fvar<var>>>::value;
  EXPECT_FALSE(val);

  val = is_stan_scalar_container<std::vector<double>>::value;
  EXPECT_TRUE(val);
  val = is_stan_scalar_container<std::vector<var>>::value;
  EXPECT_TRUE(val);
  val = is_stan_scalar_container<std::vector<fvar<double>>>::value;
  EXPECT_FALSE(val);
  val = is_stan_scalar_container<std::vector<fvar<var>>>::value;
  EXPECT_FALSE(val);
  val = is_stan_scalar_container<std::vector<fvar<fvar<double>>>>::value;
  EXPECT_FALSE(val);
  val = is_stan_scalar_container<std::vector<fvar<fvar<var>>>>::value;
  EXPECT_FALSE(val);

  val = is_stan_scalar_container<Eigen::Matrix<double, -1, -1>>::value;
  EXPECT_TRUE(val);
  val = is_stan_scalar_container<Eigen::Matrix<var, -1, -1>>::value;
  EXPECT_TRUE(val);
  val = is_stan_scalar_container<Eigen::Matrix<fvar<double>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_stan_scalar_container<Eigen::Matrix<fvar<var>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_stan_scalar_container<
      Eigen::Matrix<fvar<fvar<double>>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_stan_scalar_container<Eigen::Matrix<fvar<fvar<var>>, -1, -1>>::value;
  EXPECT_FALSE(val);
}

TEST(is_check, is_contains_ad_type) {
  using stan::is_ad_type_container;
  using stan::math::var;

  auto val = is_ad_type_container<double>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<var>::value;
  EXPECT_TRUE(val);
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
  EXPECT_TRUE(val);
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
  EXPECT_TRUE(val);
  val = is_ad_type_container<Eigen::Matrix<fvar<double>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<Eigen::Matrix<fvar<var>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<Eigen::Matrix<fvar<fvar<double>>, -1, -1>>::value;
  EXPECT_FALSE(val);
  val = is_ad_type_container<Eigen::Matrix<fvar<fvar<var>>, -1, -1>>::value;
  EXPECT_FALSE(val);
}

TEST(is_check, is_eigen_or_stan_scalar) {
  using stan::is_eigen_or_stan_scalar;
  using stan::math::var;
  auto val = is_eigen_or_stan_scalar<double>::value;
  EXPECT_TRUE(val);
  val = is_eigen_or_stan_scalar<var>::value;
  EXPECT_TRUE(val);
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

TEST(is_check, is_eigen_arithmetic) {
  using stan::is_eigen_arithmetic;
  using stan::math::var;

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
  using stan::math::var;

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
  EXPECT_TRUE(val);
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
  using stan::math::var;

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
