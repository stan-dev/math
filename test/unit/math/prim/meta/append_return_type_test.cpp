#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(MathMetaPrim, test_append_return_type) {
  test::expect_same_type<int, stan::math::append_return_type<int, int>::type>();
  test::expect_same_type<
      double, stan::math::append_return_type<double, double>::type>();
  test::expect_same_type<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
      stan::math::append_return_type<
          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
          Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>::type>();
}
