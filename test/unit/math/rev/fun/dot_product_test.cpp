#include <stan/math/rev/core.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/math/rev/test_var_value_helper.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

TEST(MathRev, scalar_vector) {
  Eigen::VectorXd vec = Eigen::VectorXd::Random(3);
  Eigen::RowVectorXd rvec = Eigen::RowVectorXd::Random(3);

  auto f1 = [&](const auto& arg) {
    return stan::math::dot_product(vec, arg) + stan::math::dot_product(arg, vec);
  };

  auto f2 = [&](const auto& arg) {
    return stan::math::dot_product(rvec, arg) + stan::math::dot_product(arg, rvec);
  };

  expect_ad2(f1, vec);
  expect_ad2(f2, rvec);
}
