#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, arena_matrix_matrix_test) {
  using Eigen::MatrixXd;
  using stan::math::arena_matrix;
  auto&& var_alloc_stack_ref = stan::math::ChainableStack::instance_->var_alloc_stack_;
  // construction
  arena_matrix<MatrixXd> a;
  arena_matrix<MatrixXd> a2;
  arena_matrix<MatrixXd> a3;
  arena_matrix<MatrixXd> b(3, 2);
  arena_matrix<MatrixXd> b2(4, 5);
  arena_matrix<MatrixXd> c(MatrixXd::Ones(3, 2));
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));
  arena_matrix<MatrixXd> d(c);
  arena_matrix<MatrixXd> e(2 * d);
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));

  // assignment
  a = c;
  a2 = std::move(d);
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));
  a3 = 2 * a;
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));
  b = d;
  b2 = std::move(c);
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));
  e = e + a;
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));
  a = MatrixXd::Ones(3, 2);
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));
  Eigen::Matrix<stan::math::var, -1, -1> plain_matvar = Eigen::MatrixXd::Ones(3, 2);
  using stan::math::add;
  // 1 move + 5 calls returning back rvalues that will be moved
  Eigen::Matrix<stan::math::var, -1, -1> ret = add(add(add(add(add(add(std::move(plain_matvar), a), a2), a3), b), b2), e);
  EXPECT_MATRIX_EQ(ret.val(), Eigen::MatrixXd::Ones(3, 2) * 10);
  ret.sum().grad();
  EXPECT_EQ(6, (var_alloc_stack_ref.size()));

  stan::math::recover_memory();
}

TEST(AgradRev, arena_matrix_vector_test) {
  using Eigen::VectorXd;
  using stan::math::arena_matrix;
  auto&& var_alloc_stack_ref = stan::math::ChainableStack::instance_->var_alloc_stack_;

  // construction
  arena_matrix<VectorXd> a;
  arena_matrix<VectorXd> a2;
  arena_matrix<VectorXd> a3;
  arena_matrix<VectorXd> b(3);
  arena_matrix<VectorXd> b2(4);
  arena_matrix<VectorXd> c(VectorXd::Ones(3));
  arena_matrix<VectorXd> d(c);
  arena_matrix<VectorXd> e(2 * d);
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));

  // assignment
  a = c;
  a2 = std::move(d);
  a3 = 2 * a;
  b = d;
  b2 = std::move(c);
  e = e + a;
  a = VectorXd::Ones(3);
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));

  Eigen::Matrix<stan::math::var, -1, 1> plain_matvar = Eigen::VectorXd::Ones(3);
  using stan::math::add;
  // 1 move + 5 calls returning back rvalues that will be moved
  Eigen::Matrix<stan::math::var, -1, 1> ret = add(add(add(add(add(add(std::move(plain_matvar), a), a2), a3), b), b2), e);
  EXPECT_MATRIX_EQ(ret.val(), Eigen::VectorXd::Ones(3) * 10);
  ret.sum().grad();
  EXPECT_EQ(6, (var_alloc_stack_ref.size()));
  stan::math::recover_memory();
}

TEST(AgradRev, arena_matrix_row_vector_test) {
  using Eigen::RowVectorXd;
  using stan::math::arena_matrix;
  auto&& var_alloc_stack_ref = stan::math::ChainableStack::instance_->var_alloc_stack_;

  // construction
  arena_matrix<RowVectorXd> a;
  arena_matrix<RowVectorXd> a2;
  arena_matrix<RowVectorXd> a3;
  arena_matrix<RowVectorXd> b(3);
  arena_matrix<RowVectorXd> b2(4);
  arena_matrix<RowVectorXd> c(RowVectorXd::Ones(3));
  arena_matrix<RowVectorXd> d(c);
  arena_matrix<RowVectorXd> e(2 * d);
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));

  // assignment
  a = c;
  a2 = std::move(d);
  a3 = 2 * a;
  b = d;
  b2 = std::move(c);
  e = e + a;
  a = RowVectorXd::Ones(3);
  EXPECT_EQ(0, (var_alloc_stack_ref.size()));

  Eigen::Matrix<stan::math::var, 1, -1> plain_matvar = Eigen::VectorXd::Ones(3);
  using stan::math::add;
  // 1 move + 5 calls returning back rvalues that will be moved
  Eigen::Matrix<stan::math::var, 1, -1> ret = add(add(add(add(add(add(std::move(plain_matvar), a), a2), a3), b), b2), e);
  EXPECT_MATRIX_EQ(ret.val(), Eigen::RowVectorXd::Ones(3) * 10);
  ret.sum().grad();
  EXPECT_EQ(6, (var_alloc_stack_ref.size()));
  stan::math::recover_memory();
}

TEST(AgradRev, arena_matrix_transpose_test) {
  using stan::math::arena_matrix;

  Eigen::VectorXd c = Eigen::VectorXd::Random(3);
  Eigen::RowVectorXd d = Eigen::RowVectorXd::Random(3);

  // The VectorXd/RowVectorXd mixup here is on purpose to test a vector can
  // initialize a row vector and visa versa
  arena_matrix<Eigen::VectorXd> a = d;
  arena_matrix<Eigen::RowVectorXd> b = c;
  EXPECT_MATRIX_EQ(a, d.transpose());
  EXPECT_MATRIX_EQ(b, c.transpose());

  EXPECT_NO_THROW(a = b);
  EXPECT_MATRIX_EQ(a, b.transpose());
  a = Eigen::VectorXd::Random(3);
  EXPECT_NO_THROW(b = a);
  EXPECT_MATRIX_EQ(a, b.transpose());

  stan::math::recover_memory();
}


TEST(AgradRev, arena_matrix_rvalue_test) {
  using stan::arena_t;
  constexpr size_t initial_bytes = stan::math::internal::DEFAULT_INITIAL_NBYTES;
  // Just below the size that would force a new allocation.
  constexpr size_t initial_mat_size = initial_bytes/sizeof(stan::math::vari);
  using vec_v = Eigen::Matrix<stan::math::var, -1, 1>;
  vec_v A = Eigen::VectorXd::Random(initial_mat_size);
  EXPECT_EQ(initial_bytes, (stan::math::ChainableStack::instance_->memalloc_.bytes_allocated()));
  // Moving causes no new alloc on stack
  arena_t<vec_v> A_arena = std::move(A);
  auto&& var_alloc_stack_ref = stan::math::ChainableStack::instance_->var_alloc_stack_;
  EXPECT_EQ(initial_bytes, (stan::math::ChainableStack::instance_->memalloc_.bytes_allocated()));
  EXPECT_EQ(1, (var_alloc_stack_ref.size()));
  // Checking transpose assign move works
  Eigen::RowVectorXd C = Eigen::RowVectorXd::Random(initial_mat_size);
  arena_t<Eigen::VectorXd> C_arena = std::move(C);
  EXPECT_EQ(2, (var_alloc_stack_ref.size()));
  EXPECT_EQ(initial_bytes, (stan::math::ChainableStack::instance_->memalloc_.bytes_allocated()));
  Eigen::VectorXd B = Eigen::VectorXd::Random(initial_mat_size);
  // Scalar type mismatch so copy is performed
  arena_t<vec_v> B_arena = std::move(B);
  EXPECT_EQ(2, (var_alloc_stack_ref.size()));
  EXPECT_EQ(196608, (stan::math::ChainableStack::instance_->memalloc_.bytes_allocated()));
  stan::math::recover_memory();
}
