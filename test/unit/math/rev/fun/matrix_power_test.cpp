#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(AgradRev, matrix_power_rev) {
  using stan::math::var;

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x
      = Eigen::MatrixXd::Random(1, 1);

  stan::math::set_zero_all_adjoints();
  auto y = stan::math::matrix_power(x, 2);
  y(0, 0).grad();

  auto adj = x.adj().eval();

  stan::math::set_zero_all_adjoints();
  auto y_ref = stan::math::pow(x(0, 0), 2);
  y_ref.grad();

  EXPECT_FLOAT_EQ(y(0, 0).val(), y_ref.val());
  EXPECT_FLOAT_EQ(adj(0, 0), x.adj()(0, 0));
}
