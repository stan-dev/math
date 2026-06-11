// Regression test for Eigen 5's eager `.select()` evaluation on autodiff
// scalars. Eigen 3.4's Select expression evaluated only the chosen branch
// per coefficient; Eigen 5 reimplements select as a CwiseTernaryOp whose
// evaluator computes both branches eagerly. For taped AD types evaluation
// has a side effect (vari allocation), and an orphan vari whose chain rule
// multiplies the adjoint by a value-derived quantity (e.g. exp: adj *
// exp(val) with val = +INF) injects 0 * INF = NaN into the adjoint of a
// live input. Stan restores lazy semantics via the ternary_evaluator
// specialization in stan/math/rev/core/Eigen_NumTraits.hpp; this test
// guards that specialization (it fails with NaN adjoints without it).
#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST_F(AgradRev, eigen_select_does_not_evaluate_unselected_branch) {
  using stan::math::var;
  Eigen::Array<var, Eigen::Dynamic, 1> x(2);
  x << 1.0, std::numeric_limits<double>::infinity();

  // Guard idiom used throughout prim/prob: only exponentiate where x < 2.
  Eigen::Array<var, Eigen::Dynamic, 1> y = (x < 2.0).select(x.exp(), 0.0);

  EXPECT_DOUBLE_EQ(y(0).val(), std::exp(1.0));
  EXPECT_DOUBLE_EQ(y(1).val(), 0.0);

  var total = stan::math::sum(y);
  total.grad();

  EXPECT_DOUBLE_EQ(x(0).adj(), std::exp(1.0));
  // Lazy select: exp(x(1)) never evaluated -> adj == 0.
  // Eager select: orphan exp(INF) vari -> 0 * INF = NaN.
  EXPECT_TRUE(std::isfinite(x(1).adj())) << "x(1).adj() = " << x(1).adj();
}
