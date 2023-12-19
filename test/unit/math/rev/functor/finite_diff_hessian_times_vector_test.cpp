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

struct fun2 {
  // fun2(x, y) = (x^2 * y) + (3 * y^2) + (5 * x * y) + sin(x)
  // d/dx fun2(x, y) = (2 * x * y) + (5 * y) + cos(x)
  // d/dy fun2(x, y) = (x^2) + (6 * y) + (5 * x)
  // d^2/dx^2 fun2(x, y) = (2 * y) - (sin(x))
  // d^2/dydx fun2(x, y) = (2 * x) + 5
  // d^2/dxdy fun2(x, y) = (2 * x) + 5
  // d^2/dy^2 fun2(x, y) = 6
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    using std::sin;
    return x(0) * x(0) * x(1) + 3.0 * x(1) * x(1) + 5.0 * x(0) * x(1)
           + sin(x(0));
  }
};

TEST(RevFunctor, finiteDiffHessianTimesVector) {
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

TEST(RevFunctor, finiteDiffHessianTimesVector2) {
  using stan::math::internal::finite_diff_hessian_times_vector_auto;

  fun2 f;

  Matrix<double, Dynamic, 1> x(2);
  x << 13, -4;

  Matrix<double, Dynamic, 1> v(2);
  v << 10, 0.2;

  Matrix<double, Dynamic, 1> Hv;
  double fx;
  finite_diff_hessian_times_vector_auto(f, x, v, fx, Hv);

  EXPECT_FLOAT_EQ(13 * 13 * -4 + 3 * -4 * -4 + 5 * 13 * -4 + std::sin(13), fx);

  EXPECT_EQ(2, Hv.size());
  EXPECT_FLOAT_EQ((2 * x(1) - std::sin(x(0))) * v(0) + (2 * x(0) + 5) * v(1),
                  Hv(0));
  EXPECT_FLOAT_EQ((2 * x(0) + 5) * v(0) + 6 * v(1), Hv(1));
}

}  // namespace finite_diff_hessian_times_vector_test
