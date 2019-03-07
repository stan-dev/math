#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ErrorHandling, isMatchingSize) {
  std::vector<double> a;
  std::vector<double> b;
  EXPECT_TRUE(stan::math::is_matching_size(a, b));

  a = {3, 3};
  EXPECT_FALSE(stan::math::is_matching_size(a, b));

  b = {3, 3};
  EXPECT_TRUE(stan::math::is_matching_size(a, b));
}

TEST(ErrorHandling, isMatchingSize_arr) {
  std::array<double, 4> a;
  std::array<double, 4> b;
  EXPECT_TRUE(stan::math::is_matching_size(a, b));

  std::array<double, 0> c;
  EXPECT_FALSE(stan::math::is_matching_size(c, b));

  std::array<double, 0> d;
  EXPECT_TRUE(stan::math::is_matching_size(c, d));
}

TEST(ErrorHandling, isMatchingSize_mat) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b;
  a.resize(4, 4);
  b.resize(4, 4);
  EXPECT_TRUE(stan::math::is_matching_size(a, b));

  a.resize(2, 4);
  EXPECT_FALSE(stan::math::is_matching_size(a, b));

  b.resize(2, 4);
  EXPECT_TRUE(stan::math::is_matching_size(a, b));
}
