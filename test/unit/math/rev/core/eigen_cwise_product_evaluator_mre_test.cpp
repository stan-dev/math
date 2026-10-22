#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

// Minimal instantiation of Eigen's binary_evaluator for a coefficient-wise
// product of two var arrays. `a * a` is a
// `CwiseBinaryOp<scalar_product_op<var, var>, ...>`, and `.sum()` builds its
// evaluator directly.
TEST_F(AgradRev, eigen_cwise_product_evaluator_var_var) {
  using stan::math::var;
  Eigen::Array<var, -1, -1> a(2, 2);
  a << 1, 2, 3, 4;
  const var var_var = (a * a).sum();
  EXPECT_DOUBLE_EQ(var_var.val(), 30.0);
}
