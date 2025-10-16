#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

TEST_F(AgradRev, MetaTraitsRevScal_is_var) {
  using stan::is_var;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_TRUE(is_var<stan::math::var>::value);
  EXPECT_TRUE((is_var<stan::math::var_value<float>>::value));
  EXPECT_TRUE((is_var<stan::math::var_value<long double>>::value));
  EXPECT_FALSE(is_var<stan::math::vari>::value);
  EXPECT_FALSE((is_var<double>::value));
  EXPECT_FALSE((is_var<stan::math::vari_value<double>>::value));
}

TEST_F(AgradRev, MetaTraitsRevScal_is_any_var_scalar) {
  using stan::is_any_var_scalar;
  using stan::is_any_var_scalar_v;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_TRUE(is_any_var_scalar<stan::math::var>::value);
  EXPECT_TRUE((is_any_var_scalar<stan::math::var_value<float>>::value));
  EXPECT_TRUE((is_any_var_scalar<stan::math::var_value<long double>>::value));
  EXPECT_FALSE(is_any_var_scalar<stan::math::vari>::value);
  EXPECT_FALSE((is_any_var_scalar<double>::value));
  EXPECT_FALSE((is_any_var_scalar<stan::math::vari_value<double>>::value));

  EXPECT_TRUE((is_any_var_scalar_v<
               std::tuple<stan::math::var, double, std::vector<double>>>));
  EXPECT_TRUE((
      is_any_var_scalar_v<std::tuple<
          std::vector<Eigen::Matrix<double, -1, 1, 0, -1, 1>,
                      std::allocator<Eigen::Matrix<double, -1, 1, 0, -1, 1>>> &,
          const double &, stan::math::var_value<double, void> &>>));
  EXPECT_TRUE(
      (is_any_var_scalar_v<
          Eigen::Matrix<double, -1, 1, 0, -1, 1>,
          Eigen::Matrix<double, 0, 0, 0, 0, 0>,
          std::tuple<
              std::tuple<std::vector<Eigen::Matrix<double, -1, 1, 0, -1, 1>> &,
                         stan::math::var_value<double, void> &>,
              const double &>>));
}
