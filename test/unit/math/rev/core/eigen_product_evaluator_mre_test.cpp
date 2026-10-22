#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

// Minimal instantiation of Eigen's product_evaluator for a product of two
// var matrices. `(A * A).sum()` builds the product's evaluator directly
// (unlike `C = A * A`, which goes through Eigen's Assignment).
TEST_F(AgradRev, eigen_product_evaluator_var_var) {
  using stan::math::var;
  Eigen::Matrix<var, -1, -1> a(2, 2);
  a << 1, 2, 3, 4;
  const var var_var = (a * a).sum();
  EXPECT_DOUBLE_EQ(var_var.val(), (a.val() * a.val()).sum());
}
