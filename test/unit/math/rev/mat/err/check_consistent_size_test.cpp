#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

TEST(AgradRevErrorHandlingScalar, checkConsistentSize) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::check_consistent_size;
  using stan::size_of;
  using stan::math::var;

  const char* function = "check_consistent_size";
  const char* name1 = "name1";
  

  Matrix<var,Dynamic,1> v1(4);
  v1 << 4.0,5.0,6.0,7.0;
  EXPECT_EQ(4U, size_of(v1));
  EXPECT_NO_THROW(check_consistent_size(function, name1, v1, 4U));
  EXPECT_THROW(check_consistent_size(function, name1, v1, 2U), std::invalid_argument);
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, checkConsistentSize_nan) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::check_consistent_size;
  using stan::size_of;
  using stan::math::var;

  const char* function = "check_consistent_size";
  const char* name1 = "name1";
  
  double nan = std::numeric_limits<double>::quiet_NaN();

  Matrix<var,Dynamic,1> v1(4);
  v1 << nan,nan,4,nan;
  EXPECT_EQ(4U, size_of(v1));
  EXPECT_NO_THROW(check_consistent_size(function, name1, v1, 4U));
  EXPECT_THROW(check_consistent_size(function, name1, v1, 2U), std::invalid_argument);
  stan::math::recover_memory();
}
