#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::var;

TEST(AgradRevErrorHandlingArray, CheckPositive) {
  using stan::math::check_positive;
  const char* function = "check_positive";

  std::vector<var> x;
  x.push_back(var(1.0));
  x.push_back(var(2.0));
  x.push_back(var(3.0));

  EXPECT_NO_THROW(check_positive(function, "x", x));

  stan::math::recover_memory();
}
#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

using stan::math::var;

TEST(AgradRevErrorHandlingMatrix, CheckPositive) {
  using stan::math::check_positive;
  const char* function = "check_positive";

  Eigen::Matrix<var, Eigen::Dynamic, 1> x_mat(3);
  x_mat << 1, 2, 3;
  for (int i = 0; i < x_mat.size(); i++) {
    EXPECT_NO_THROW(check_positive(function, "x", x_mat));
  }

  x_mat(0) = 0;

  EXPECT_THROW(check_positive(function, "x", x_mat), std::domain_error);

  stan::math::recover_memory();
}
