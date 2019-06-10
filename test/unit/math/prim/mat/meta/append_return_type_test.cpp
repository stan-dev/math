#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

using stan::math::append_return_type;

TEST(MetaTraits, test_append_return_type) {
  test::expect_same_type<int, append_return_type<int, int>::type>();
  test::expect_same_type<double, append_return_type<double, double>::type>();
  test::expect_same_type<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
      append_return_type<
          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
          Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>::type>();
}
