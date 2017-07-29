#include <stan/math/fwd/mat.hpp>
#include <test/unit/math/mix/mat/fun/append_array_test.hpp>
#include <gtest/gtest.h>

TEST(MathFunctionsFwd, append_array) {
  check<int, fd, fd>();
  check<fd, int, fd>();
  check<double, fd, fd>();
  check<fd, double, fd>();
  check<fd, fd, fd>();
  check<V, Vfd, Vfd>();
  check<Vfd, V, Vfd>();
  check<Vfd, Vfd, Vfd>();
  check<RV, RVfd, RVfd>();
  check<RVfd, RV, RVfd>();
  check<RVfd, RVfd, RVfd>();
  check<M, Mfd, Mfd>();
  check<Mfd, M, Mfd>();
  check<Mfd, Mfd, Mfd>();

  check<int, ffd, ffd>();
  check<ffd, int, ffd>();
  check<double, ffd, ffd>();
  check<ffd, double, ffd>();
  check<ffd, ffd, ffd>();
  check<V, Vffd, Vffd>();
  check<Vffd, V, Vffd>();
  check<Vffd, Vffd, Vffd>();
  check<RV, RVffd, RVffd>();
  check<RVffd, RV, RVffd>();
  check<RVffd, RVffd, RVffd>();
  check<M, Mffd, Mffd>();
  check<Mffd, M, Mffd>();
  check<Mffd, Mffd, Mffd>();
}
