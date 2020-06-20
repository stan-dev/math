#include <stan/math/rev/core.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/math/rev/test_var_value_helper.hpp>
#include <stan/math/rev/fun/rows_dot_product.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

TEST(MathRev, var_value_tests) {
  Eigen::MatrixXd mat1 = Eigen::MatrixXd::Random(3, 3);
  Eigen::MatrixXd mat2 = Eigen::MatrixXd::Random(3, 3);

  auto f = [&](const auto& arg) {
    return (stan::math::rows_dot_product(mat1, arg) +
    stan::math::rows_dot_product(arg, mat1)).eval();
  };

  expect_ad2(f, mat2);
}
