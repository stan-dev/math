#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

// Minimal instantiation of Eigen's Assignment of a product of two var
// matrices. `c = a * a` assigns a `Product<...>` to a var matrix with
// `assign_op<var, var>`, which Eigen handles in `Assignment` without a
// product_evaluator.
TEST_F(AgradRev, eigen_product_assignment_var_var) {
  using stan::math::var;
  Eigen::Matrix<var, -1, -1> a(2, 2);
  a << 1, 2, 3, 4;
  Eigen::Matrix<var, -1, -1> c(2, 2);
  c = a * a;
  EXPECT_DOUBLE_EQ(c.val().sum(), 54.0);
}
