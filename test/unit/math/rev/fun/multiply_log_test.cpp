#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(RevTest, multiply_log_issue3146) {
  std::int32_t x = 1;
  using stan::math::var_value;
  using mat_t = Eigen::Matrix<double, -1, -1, 1, -1, -1>;
  var_value<mat_t> y = Eigen::MatrixXd::Random(10, 10);
  auto z = stan::math::multiply_log(x, y);
}
