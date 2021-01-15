#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>

#include <gtest/gtest.h>

TEST(MathFunRev, increment_adjoint_non_var) {
  double a = 0;
  double a1 = a;
  Eigen::VectorXd b(3);
  b << 1, 2, 3;
  Eigen::VectorXd b1 = b;
  stan::math::increment_adjoint(a, 1);
  stan::math::increment_adjoint(a, b);
  stan::math::increment_adjoint(b, 2);
  stan::math::increment_adjoint(b, b);
  EXPECT_EQ(a, a1);
  EXPECT_MATRIX_EQ(b, b1);
}

TEST(MathFunRev, increment_adjoint_scal_scal) {
  stan::math::var a = 0;
  a.adj() = 2;
  stan::math::increment_adjoint(a, 1);
  EXPECT_EQ(a.adj(), 3);
}

TEST(MathFunRev, increment_adjoint_scal_vec) {
  stan::math::var a = 0;
  a.adj() = 2;
  Eigen::VectorXd b(3);
  b << 1, 2, 3;
  stan::math::increment_adjoint(a, b);
  EXPECT_EQ(a.adj(), 8);
}

TEST(MathFunRev, increment_adjoint_vec_scal) {
  Eigen::VectorXd adj1(3);
  adj1 << 4, 5, 6;
  Eigen::Matrix<stan::math::var, -1, 1> a(3);
  a << 1, 2, 3;
  a.adj() = adj1;

  stan::math::increment_adjoint(a, 3);
  EXPECT_MATRIX_EQ(a.adj(), adj1.array() + 3);
}

TEST(MathFunRev, increment_adjoint_vec_vec) {
  Eigen::VectorXd adj1(3);
  adj1 << 4, 5, 6;
  Eigen::VectorXd adj2(3);
  adj2 << 7, 8, 9;
  Eigen::Matrix<stan::math::var, -1, 1> a(3);
  a << 1, 2, 3;
  a.adj() = adj1;

  stan::math::increment_adjoint(a, adj2);
  EXPECT_MATRIX_EQ(a.adj(), adj1 + adj2);
}
