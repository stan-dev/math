#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(MathMetaPrim, test_append_return_type) {
  EXPECT_SAME_TYPE(int, stan::math::append_return_type<int, int>::type);
  EXPECT_SAME_TYPE(double,
                   stan::math::append_return_type<double, double>::type);
  EXPECT_SAME_TYPE(
      Eigen::MatrixXd,
      stan::math::append_return_type<Eigen::MatrixXd, Eigen::MatrixXi>::type);
}
