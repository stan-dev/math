#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, arena_matrix_matrix_test) {
  using Eigen::MatrixXd;
  using stan::math::arena_matrix;

  // construction
  arena_matrix<MatrixXd> a;
  arena_matrix<MatrixXd> a2;
  arena_matrix<MatrixXd> a3;
  arena_matrix<MatrixXd> b(3, 2);
  arena_matrix<MatrixXd> b2(4, 5);
  arena_matrix<MatrixXd> c(MatrixXd::Ones(3, 2));
  arena_matrix<MatrixXd> d(c);
  arena_matrix<MatrixXd> e(2 * d);

  // assignment
  a = c;
  a2 = std::move(d);
  a3 = 2 * a;
  b = d;
  b2 = std::move(c);
  e = e + a;
  a = MatrixXd::Ones(3, 2);

  EXPECT_MATRIX_EQ(a + a2 + a3 + b + b2 + e, MatrixXd::Ones(3, 2) * 9);
  stan::math::recover_memory();
}

TEST(AgradRev, arena_matrix_vector_test) {
  using Eigen::VectorXd;
  using stan::math::arena_matrix;

  // construction
  arena_matrix<VectorXd> a;
  arena_matrix<VectorXd> a2;
  arena_matrix<VectorXd> a3;
  arena_matrix<VectorXd> b(3);
  arena_matrix<VectorXd> b2(4);
  arena_matrix<VectorXd> c(VectorXd::Ones(3));
  arena_matrix<VectorXd> d(c);
  arena_matrix<VectorXd> e(2 * d);

  // assignment
  a = c;
  a2 = std::move(d);
  a3 = 2 * a;
  b = d;
  b2 = std::move(c);
  e = e + a;
  a = VectorXd::Ones(3);

  EXPECT_MATRIX_EQ(a + a2 + a3 + b + b2 + e, VectorXd::Ones(3) * 9);
  stan::math::recover_memory();
}

TEST(AgradRev, arena_matrix_row_vector_test) {
  using Eigen::RowVectorXd;
  using stan::math::arena_matrix;

  // construction
  arena_matrix<RowVectorXd> a;
  arena_matrix<RowVectorXd> a2;
  arena_matrix<RowVectorXd> a3;
  arena_matrix<RowVectorXd> b(3);
  arena_matrix<RowVectorXd> b2(4);
  arena_matrix<RowVectorXd> c(RowVectorXd::Ones(3));
  arena_matrix<RowVectorXd> d(c);
  arena_matrix<RowVectorXd> e(2 * d);

  // assignment
  a = c;
  a2 = std::move(d);
  a3 = 2 * a;
  b = d;
  b2 = std::move(c);
  e = e + a;
  a = RowVectorXd::Ones(3);

  EXPECT_MATRIX_EQ(a + a2 + a3 + b + b2 + e, RowVectorXd::Ones(3) * 9);
  stan::math::recover_memory();
}

TEST(AgradRev, arena_matrix_auto_transpose) {
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::arena_matrix;

  VectorXd v1(2);
  v1 << 1, 2;
  VectorXd v2(2);
  v2 << 3, 4;
  RowVectorXd r1 = v1;
  RowVectorXd r2 = v2;

  arena_matrix<RowVectorXd> ra(v1);
  EXPECT_MATRIX_EQ(ra, r1);
  ra = v2;
  EXPECT_MATRIX_EQ(ra, r2);

  arena_matrix<VectorXd> va(r1);
  EXPECT_MATRIX_EQ(va, v1);
  va = r2;
  EXPECT_MATRIX_EQ(va, v2);
}
