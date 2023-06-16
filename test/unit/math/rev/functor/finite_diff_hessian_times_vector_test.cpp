#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <iostream>
#include <stdexcept>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;

namespace finite_diff_hessian_times_vector_test {
// fun1(x, y) = (x^2 * y) + (3 * y^2)
struct fun1 {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return x(0) * x(0) * x(1) + 3.0 * x(1) * x(1);
  }
};

TEST(MixFunctor, hessianTimesVector) {
  using stan::math::internal::finite_diff_hessian_times_vector_auto;

  fun1 f;

  Matrix<double, Dynamic, 1> x(2);
  x << 2, -3;

  Matrix<double, Dynamic, 1> v(2);
  v << 8, 5;

  Matrix<double, Dynamic, 1> Hv;
  double fx;
  finite_diff_hessian_times_vector_auto(f, x, v, fx, Hv);

  EXPECT_FLOAT_EQ(2 * 2 * -3 + 3.0 * -3 * -3, fx);

  EXPECT_EQ(2, Hv.size());
  EXPECT_FLOAT_EQ(2 * x(1) * v(0) + 2 * x(0) * v(1), Hv(0));
  EXPECT_FLOAT_EQ(2 * x(0) * v(0) + 6 * v(1), Hv(1));
}
}  // namespace finite_diff_hessian_times_vector_test
