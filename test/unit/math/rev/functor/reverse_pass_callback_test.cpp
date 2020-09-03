#include <stan/math.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, reverse_pass_callback_test) {
  stan::math::var a = 1;
  stan::math::var b = 1;

  stan::math::reverse_pass_callback(
      [=]() mutable { a.vi_->adj_ = b.adj() * 3; });
  EXPECT_EQ(a.adj(), 0);
  b.grad();
  EXPECT_EQ(a.adj(), 3);
}
