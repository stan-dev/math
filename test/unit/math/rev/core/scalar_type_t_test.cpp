#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>

TEST(AgradRev, scalar_type_t_var_value) {
  EXPECT_TRUE(
      (std::is_same<
          stan::math::var,
          stan::scalar_type_t<stan::math::var_value<Eigen::MatrixXd>>>::value));

  EXPECT_TRUE((
      std::is_same<stan::scalar_type_t<stan::math::var_value<Eigen::MatrixXd>>,
                   stan::scalar_type_t<
                       const stan::math::var_value<Eigen::MatrixXd>&>>::value));
}
