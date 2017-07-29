#include <stan/math/rev/mat.hpp>
#include <test/unit/math/mix/mat/fun/append_array_test.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, append_array_rev) {
  check<int, var, var>();
  check<var, int, var>();
  check<double, var, var>();
  check<var, double, var>();
  check<var, var, var>();
  check<V, Vv, Vv>();
  check<Vv, V, Vv>();
  check<Vv, Vv, Vv>();
  check<RV, RVv, RVv>();
  check<RVv, RV, RVv>();
  check<RVv, RVv, RVv>();
  check<M, Mv, Mv>();
  check<Mv, M, Mv>();
  check<Mv, Mv, Mv>();
}
