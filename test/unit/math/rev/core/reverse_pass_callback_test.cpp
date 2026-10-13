#include <stan/math.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

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
  // Fixed-size Eigen matrices require more than 8 byte alignment, and a
  // callback stores its captures inside the vari allocated on the arena.
  const Eigen::Matrix2d m = Eigen::Matrix2d::Identity();
  stan::math::var b = 1;
  bool captured_aligned = false;
  double captured_sum = 0;
  // leave the next free arena byte 8 but not 32 byte aligned
  stan::math::ChainableStack::instance_->memalloc_.alloc(8);
  stan::math::reverse_pass_callback([m, &captured_aligned, &captured_sum]() {
    captured_aligned
        = stan::math::is_aligned(m.data(), alignof(Eigen::Matrix2d));
    captured_sum = (m * m).sum();
  });
  b.grad();
  EXPECT_TRUE(captured_aligned);
  EXPECT_FLOAT_EQ(captured_sum, 2.0);
}
