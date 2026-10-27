#include <stan/math.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

namespace {
/**
 * Leave the next free byte of the arena 8 bytes past a 16 byte boundary.
 */
void misalign_arena() {
  auto& arena = stan::math::ChainableStack::instance_->memalloc_;
  while (!stan::math::is_aligned(arena.alloc(8), 16U)) {
  }
}
}  // namespace

TEST_F(AgradRev, Rev_reverse_pass_callback_test) {
  stan::math::var a = 1;
  stan::math::var b = 1;

  stan::math::reverse_pass_callback(
      [=]() mutable { a.vi_->adj_ = b.adj() * 3; });
  EXPECT_EQ(a.adj(), 0);
  b.grad();
  EXPECT_EQ(a.adj(), 3);
}

TEST_F(AgradRev, Rev_reverse_pass_callback_over_aligned_capture) {
  // a callback stores its captures inside the vari allocated on the arena
  const Eigen::Vector2d m(1.0, 2.0);
  stan::math::var b = 1;
  bool captured_aligned = false;
  double captured_sum = 0;
  misalign_arena();
  stan::math::reverse_pass_callback([m, &captured_aligned, &captured_sum]() {
    captured_aligned
        = stan::math::is_aligned(m.data(), alignof(Eigen::Vector2d));
    captured_sum = m.sum();
  });
  b.grad();
  EXPECT_TRUE(captured_aligned);
  EXPECT_FLOAT_EQ(captured_sum, 3.0);
}

TEST_F(AgradRev, Rev_multiply_fixed_size_var_matrices) {
  using stan::math::var;
  Eigen::Matrix<var, 2, 2> f;
  f << 1, 2, 3, 4;
  misalign_arena();
  Eigen::Matrix<var, 2, 2> c = stan::math::multiply(f, f);
  stan::math::sum(c).grad();
  const Eigen::Matrix2d f_val = f.val();
  EXPECT_MATRIX_NEAR(c.val(), f_val * f_val, 1e-12);
  // d sum(F * F) / dF = ones * F^T + F^T * ones
  const Eigen::Matrix2d ones = Eigen::Matrix2d::Ones();
  EXPECT_MATRIX_NEAR(
      f.adj(), ones * f_val.transpose() + f_val.transpose() * ones, 1e-12);
}
