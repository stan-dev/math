#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRevErrorHandlingScalar, checkConsistentSizes) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::check_consistent_sizes;
  using stan::math::var;
  using stan::size_of;

  const char* function = "testConsSizes";
  const char* name1 = "name1";
  const char* name2 = "name2";
  const char* name3 = "name3";
  const char* name4 = "name4";

  Matrix<var, Dynamic, 1> v1(4);
  Matrix<var, Dynamic, 1> v2(4);
  Matrix<var, Dynamic, 1> v3(4);
  Matrix<var, Dynamic, 1> v4(4);
  ASSERT_EQ(4U, size_of(v1));
  ASSERT_EQ(4U, size_of(v2));
  ASSERT_EQ(4U, size_of(v3));
  ASSERT_EQ(4U, size_of(v4));
  EXPECT_NO_THROW(check_consistent_sizes(function, name1, v1, name2, v2));
  EXPECT_NO_THROW(
      check_consistent_sizes(function, name1, v1, name2, v2, name3, v3));
  EXPECT_NO_THROW(check_consistent_sizes(function, name1, v1, name2, v2, name3,
                                         v3, name4, v4));

  Matrix<var, Dynamic, 1> v(3);

  ASSERT_EQ(3U, size_of(v));
  const char* name = "inconsistent";
  EXPECT_THROW(check_consistent_sizes(function, name, v, name2, v2),
               std::invalid_argument);
  EXPECT_THROW(check_consistent_sizes(function, name1, v1, name, v),
               std::invalid_argument);
  EXPECT_THROW(check_consistent_sizes(function, name, v, name2, v2, name3, v3),
               std::invalid_argument);
  EXPECT_THROW(check_consistent_sizes(function, name1, v1, name, v, name3, v3),
               std::invalid_argument);
  EXPECT_THROW(check_consistent_sizes(function, name1, v1, name2, v2, name, v),
               std::invalid_argument);

  EXPECT_THROW(check_consistent_sizes(function, name, v, name2, v2, name3, v3,
                                      name4, v4),
               std::invalid_argument);
  EXPECT_THROW(check_consistent_sizes(function, name1, v1, name, v, name3, v3,
                                      name4, v4),
               std::invalid_argument);
  EXPECT_THROW(check_consistent_sizes(function, name1, v1, name2, v2, name, v,
                                      name4, v4),
               std::invalid_argument);
  EXPECT_THROW(
      check_consistent_sizes(function, name, v, name2, v2, name3, v3, name, v),
      std::invalid_argument);
  stan::math::recover_memory();
}
