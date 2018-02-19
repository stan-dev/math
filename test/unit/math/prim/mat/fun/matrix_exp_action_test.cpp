#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_action.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <algorithm>
#include <random>

TEST(MathMatrix, matrix_exp_action_diag) {
  using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorType = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  stan::math::matrix_exp_action_handler handler;

  {
    double t = 0.5;
    MatrixType m1(2, 2);
    VectorType b(2);
    m1 << 2, 0, 0, 2;
    b << 1, 1;
    auto res = handler.perform_action(m1, b, t);
    EXPECT_NEAR(res(0), M_E, 1.e-8);
    EXPECT_NEAR(res(1), M_E, 1.e-8);
  }

  {
    double t = 1.0;
    MatrixType m1(2, 2);
    VectorType b = VectorType::Random(2);
    m1 << 1, 0, 0, 2;
    auto res = handler.perform_action(m1, b, t);
    EXPECT_NEAR(res(0), b(0)*M_E, 1.e-8);
    EXPECT_NEAR(res(1), b(1)*M_E*M_E, 1.e-8);
  }
}
