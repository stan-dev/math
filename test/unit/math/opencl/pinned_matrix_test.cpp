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
  pinned_matrix<MatrixXd> a2(a);
  pinned_matrix<MatrixXd> a3(std::move(a));
  pinned_matrix<MatrixXd> b(4, 5);
  pinned_matrix<MatrixXd> b2(4, 5);
  pinned_matrix<MatrixXd> b3(4, 5);
  pinned_matrix<MatrixXd> c(MatrixXd::Ones(3, 2));
  pinned_matrix<MatrixXd> d(c);
  pinned_matrix<MatrixXd> e(std::move(c));
  pinned_matrix<MatrixXd> f(2 * d);
  pinned_matrix<MatrixXd> g(0, 0);

  // assignment
  a = e;
  a2 = std::move(d);
  a3 = 2 * a;
  b = a;
  b2 = std::move(a);
  b3 = 2 * a2;
  e = e + a2;
  e += a2;

  EXPECT_MATRIX_EQ(a, MatrixXd(0, 0));
  EXPECT_MATRIX_EQ(a2, MatrixXd::Ones(3, 2));
  EXPECT_MATRIX_EQ(a3, MatrixXd::Ones(3, 2) * 2);
  EXPECT_MATRIX_EQ(b, MatrixXd::Ones(3, 2));
  EXPECT_MATRIX_EQ(b2, MatrixXd::Ones(3, 2));
  EXPECT_MATRIX_EQ(b3, MatrixXd::Ones(3, 2) * 2);
  EXPECT_MATRIX_EQ(c, MatrixXd(0, 0));
  EXPECT_MATRIX_EQ(d, MatrixXd(0, 0));
  EXPECT_MATRIX_EQ(e, MatrixXd::Ones(3, 2) * 3);
  EXPECT_MATRIX_EQ(f, MatrixXd::Ones(3, 2) * 2);

  stan::math::recover_memory();
}

TEST(AgradRev, pinned_matrix_vector_test) {
  using Eigen::VectorXd;
  using stan::math::pinned_matrix;

  // construction
  pinned_matrix<VectorXd> a;
  pinned_matrix<VectorXd> a2(a);
  pinned_matrix<VectorXd> a3(std::move(a));
  pinned_matrix<VectorXd> b(4);
  pinned_matrix<VectorXd> b2(4);
  pinned_matrix<VectorXd> b3(4);
  pinned_matrix<VectorXd> c(VectorXd::Ones(3));
  pinned_matrix<VectorXd> d(c);
  pinned_matrix<VectorXd> e(std::move(c));
  pinned_matrix<VectorXd> f(2 * d);
  pinned_matrix<VectorXd> g(0);

  // assignment
  a = e;
  a2 = std::move(d);
  a3 = 2 * a;
  b = a;
  b2 = std::move(a);
  b3 = 2 * a2;
  e = e + a2;
  e += a2;

  EXPECT_MATRIX_EQ(a, VectorXd(0));
  EXPECT_MATRIX_EQ(a2, VectorXd::Ones(3));
  EXPECT_MATRIX_EQ(a3, VectorXd::Ones(3) * 2);
  EXPECT_MATRIX_EQ(b, VectorXd::Ones(3));
  EXPECT_MATRIX_EQ(b2, VectorXd::Ones(3));
  EXPECT_MATRIX_EQ(b3, VectorXd::Ones(3) * 2);
  EXPECT_MATRIX_EQ(c, VectorXd(0));
  EXPECT_MATRIX_EQ(d, VectorXd(0));
  EXPECT_MATRIX_EQ(e, VectorXd::Ones(3) * 3);
  EXPECT_MATRIX_EQ(f, VectorXd::Ones(3) * 2);

  stan::math::recover_memory();
}

TEST(AgradRev, pinned_matrix_row_vector_test) {
  using Eigen::RowVectorXd;
  using stan::math::pinned_matrix;

  // construction
  pinned_matrix<RowVectorXd> a;
  pinned_matrix<RowVectorXd> a2(a);
  pinned_matrix<RowVectorXd> a3(std::move(a));
  pinned_matrix<RowVectorXd> b(4);
  pinned_matrix<RowVectorXd> b2(4);
  pinned_matrix<RowVectorXd> b3(4);
  pinned_matrix<RowVectorXd> c(RowVectorXd::Ones(3));
  pinned_matrix<RowVectorXd> d(c);
  pinned_matrix<RowVectorXd> e(std::move(c));
  pinned_matrix<RowVectorXd> f(2 * d);
  pinned_matrix<RowVectorXd> g(0);

  // assignment
  a = e;
  a2 = std::move(d);
  a3 = 2 * a;
  b = a;
  b2 = std::move(a);
  b3 = 2 * a2;
  e = e + a2;
  e += a2;

  EXPECT_MATRIX_EQ(a, RowVectorXd(0));
  EXPECT_MATRIX_EQ(a2, RowVectorXd::Ones(3));
  EXPECT_MATRIX_EQ(a3, RowVectorXd::Ones(3) * 2);
  EXPECT_MATRIX_EQ(b, RowVectorXd::Ones(3));
  EXPECT_MATRIX_EQ(b2, RowVectorXd::Ones(3));
  EXPECT_MATRIX_EQ(b3, RowVectorXd::Ones(3) * 2);
  EXPECT_MATRIX_EQ(c, RowVectorXd(0));
  EXPECT_MATRIX_EQ(d, RowVectorXd(0));
  EXPECT_MATRIX_EQ(e, RowVectorXd::Ones(3) * 3);

  stan::math::recover_memory();
}

#endif
