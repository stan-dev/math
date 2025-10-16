#include <stan/math/rev/core.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(AgradRev, varmat_unary_negative) {
  stan::math::var_value<Eigen::MatrixXd> x = Eigen::MatrixXd::Random(2, 3);

  stan::math::var_value<Eigen::MatrixXd> y = -x;

  EXPECT_MATRIX_EQ(y.val(), -x.val());
  for (size_t i = 0; i < y.rows(); ++i) {
    for (size_t j = 0; j < y.cols(); ++j) {
      y(i, j).grad();
      Eigen::MatrixXd adj_ref = Eigen::MatrixXd::Zero(y.rows(), y.cols());
      adj_ref(i, j) = -1.0;

      EXPECT_MATRIX_EQ(adj_ref, x.adj());
      stan::math::set_zero_all_adjoints();
    }
  }
}
