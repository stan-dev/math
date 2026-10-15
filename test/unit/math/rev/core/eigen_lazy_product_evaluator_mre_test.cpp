#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

// Minimal instantiation of Eigen's product_evaluator for a lazy product of
// two var matrices. `A.lazyProduct(A)` is a `Product<..., LazyProduct>`, and
// `.sum()` builds its evaluator directly. Stan does not specialize lazy
// products, so this is Eigen's coefficient-based evaluator.
TEST_F(AgradRev, eigen_lazy_product_evaluator_var_var) {
  using stan::math::var;
  Eigen::Matrix<var, -1, -1> a(2, 2);
  a << 1, 2, 3, 4;
  const var var_var = a.lazyProduct(a).sum();
  EXPECT_DOUBLE_EQ(var_var.val(), (a.val() * a.val()).sum());
}
