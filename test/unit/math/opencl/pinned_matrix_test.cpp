#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, pinned_matrix_matrix_test) {
  using Eigen::MatrixXd;
  using stan::math::pinned_matrix;

  // construction
  pinned_matrix<MatrixXd> a;
  pinned_matrix<MatrixXd> a2;
  pinned_matrix<MatrixXd> a3;
  pinned_matrix<MatrixXd> b(3, 2);
  pinned_matrix<MatrixXd> b2(4, 5);
  pinned_matrix<MatrixXd> c(MatrixXd::Ones(3, 2));
  pinned_matrix<MatrixXd> d(c);
  pinned_matrix<MatrixXd> e(2 * d);

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

TEST(AgradRev, pinned_matrix_vector_test) {
  using Eigen::VectorXd;
  using stan::math::pinned_matrix;

  // construction
  pinned_matrix<VectorXd> a;
  pinned_matrix<VectorXd> a2;
  pinned_matrix<VectorXd> a3;
  pinned_matrix<VectorXd> b(3);
  pinned_matrix<VectorXd> b2(4);
  pinned_matrix<VectorXd> c(VectorXd::Ones(3));
  pinned_matrix<VectorXd> d(c);
  pinned_matrix<VectorXd> e(2 * d);

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

TEST(AgradRev, pinned_matrix_row_vector_test) {
  using Eigen::RowVectorXd;
  using stan::math::pinned_matrix;

  // construction
  pinned_matrix<RowVectorXd> a;
  pinned_matrix<RowVectorXd> a2;
  pinned_matrix<RowVectorXd> a3;
  pinned_matrix<RowVectorXd> b(3);
  pinned_matrix<RowVectorXd> b2(4);
  pinned_matrix<RowVectorXd> c(RowVectorXd::Ones(3));
  pinned_matrix<RowVectorXd> d(c);
  pinned_matrix<RowVectorXd> e(2 * d);

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

#endif
