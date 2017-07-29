#include <stan/math/mix/mat.hpp>
#include <test/unit/math/mix/mat/fun/append_array_test.hpp>
#include <gtest/gtest.h>

TEST(MathFunctionsMix, append_array) {
  check<int, fv, fv>();
  check<fv, int, fv>();
  check<double, fv, fv>();
  check<fv, double, fv>();
  check<fv, fv, fv>();
  check<V, Vfv, Vfv>();
  check<Vfv, V, Vfv>();
  check<Vfv, Vfv, Vfv>();
  check<RV, RVfv, RVfv>();
  check<RVfv, RV, RVfv>();
  check<RVfv, RVfv, RVfv>();
  check<M, Mfv, Mfv>();
  check<Mfv, M, Mfv>();
  check<Mfv, Mfv, Mfv>();

  check<int, ffv, ffv>();
  check<ffv, int, ffv>();
  check<double, ffv, ffv>();
  check<ffv, double, ffv>();
  check<ffv, ffv, ffv>();
  check<V, Vffv, Vffv>();
  check<Vffv, V, Vffv>();
  check<Vffv, Vffv, Vffv>();
  check<RV, RVffv, RVffv>();
  check<RVffv, RV, RVffv>();
  check<RVffv, RVffv, RVffv>();
  check<M, Mffv, Mffv>();
  check<Mffv, M, Mffv>();
  check<Mffv, Mffv, Mffv>();
}
