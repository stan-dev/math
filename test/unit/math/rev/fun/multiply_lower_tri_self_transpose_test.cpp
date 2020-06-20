#include <stan/math/rev/core.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/math/rev/test_var_value_helper.hpp>
#include <stan/math/rev/fun/multiply_lower_tri_self_transpose.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

TEST(MathRev, var_value_tests) {
  Eigen::MatrixXd mat = Eigen::MatrixXd::Random(3, 3);

  auto f = [&](const auto& arg) {
    return stan::math::multiply_lower_tri_self_transpose(arg);
  };

  expect_ad2(f, mat);
}
