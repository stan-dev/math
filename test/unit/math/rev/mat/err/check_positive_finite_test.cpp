#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

using stan::math::check_positive_finite;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar, CheckPositiveFinite_Matrix) {
  const char* function = "check_positive_finite";
  Eigen::Matrix<var, Eigen::Dynamic, 1> x;

  x.resize(3);
  x << 3, 2, 1;
  ASSERT_NO_THROW(check_positive_finite(function, "x", x))
      << "check_positive_finite should be true with finite x";

  x.resize(3);
  x << 2, 1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on Inf";

  x.resize(3);
  x << 0, 1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on x=0";

  x.resize(3);
  x << -1, 1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on x=-1";

  x.resize(3);
  x << 2, 1, -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on -Inf";

  x.resize(3);
  x << 1, 2, std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on NaN";
  stan::math::recover_memory();
}
