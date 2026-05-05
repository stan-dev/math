// Regression test for Eigen 5's eager `.select()` evaluation on autodiff
// scalars. Eigen 3.4's Select expression evaluated only the chosen branch
// per coefficient; Eigen 5 reimplements select as a CwiseTernaryOp whose
// evaluator computes both branches eagerly. For taped AD types evaluation
// has a side effect (vari allocation), and an orphan vari whose chain rule
// multiplies the adjoint by a value-derived quantity (e.g. exp: adj *
// exp(val) with val = +INF) injects 0 * INF = NaN into the adjoint of a
// live input. Stan restores lazy semantics via the ternary_evaluator
// specializations in stan/math/{rev/core,fwd/fun}/Eigen_NumTraits.hpp.
#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathMixCore, eigen_select_does_not_evaluate_unselected_branch) {
  using stan::math::fvar;
  using stan::math::var;
  using fvar_var = fvar<var>;

  Eigen::Array<fvar_var, Eigen::Dynamic, 1> x(2);
  x(0) = fvar_var(var(1.0), 1.0);
  x(1) = fvar_var(var(std::numeric_limits<double>::infinity()), 1.0);

  // Guard idiom used throughout prim/prob: only exponentiate where x < 2.
  Eigen::Array<fvar_var, Eigen::Dynamic, 1> y = (x < 2.0).select(x.exp(), 0.0);

  EXPECT_DOUBLE_EQ(y(0).val_.val(), std::exp(1.0));
  EXPECT_DOUBLE_EQ(y(1).val_.val(), 0.0);

  fvar_var total = stan::math::sum(y);
  // d/dx[exp(x)] tangent: total.d_ == exp(x0) for the selected branch only.
  total.d_.grad();

  // Lazy select: x(1)'s exp(+INF) is never evaluated, adjoint stays 0.
  // Eager select: orphan exp(INF) vari -> 0 * INF = NaN val-side adjoint.
  EXPECT_DOUBLE_EQ(x(0).val_.adj(), std::exp(1.0));
  EXPECT_TRUE(std::isfinite(x(1).val_.adj()))
      << "x(1).val_.adj() = " << x(1).val_.adj();

  stan::math::recover_memory();
}
