#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

// as_array_or_scalar of a var_value that already holds an array returns it
// unchanged, so it can be applied twice
TEST_F(AgradRev, MathFunRev_as_array_or_scalar_var_value_array) {
  Eigen::ArrayXd a_val(3);
  a_val << 1, 2, 3;
  stan::math::var_value<Eigen::ArrayXd> a(a_val);
  auto b = stan::math::as_array_or_scalar(a);
  EXPECT_TRUE((stan::is_eigen_array<decltype(b.val())>::value));
  EXPECT_MATRIX_EQ(b.val(), a_val);

  stan::math::var_value<Eigen::VectorXd> v(a_val.matrix());
  auto c = stan::math::as_array_or_scalar(stan::math::as_array_or_scalar(v));
  EXPECT_MATRIX_EQ(c.val(), a_val);

  stan::math::var s = stan::math::make_callback_var(
      b.val().sum(), [b](auto& vi) mutable { b.adj() += vi.adj(); });
  s.grad();
  EXPECT_MATRIX_EQ(a.adj(), Eigen::ArrayXd::Ones(3));
}
