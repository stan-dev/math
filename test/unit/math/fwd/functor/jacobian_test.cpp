#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <stdexcept>

using Eigen::Dynamic;
using Eigen::Matrix;

// fun_2_3: R^2 --> R^3 | (x, y) --> [x, (x + y), (x * y)]
struct fun_2_3 {
  template <typename T>
  inline Matrix<T, Dynamic, 1> operator()(
      const Matrix<T, Dynamic, 1>& x) const {
    Matrix<T, Dynamic, 1> z(3);
    z << x(0), x(0) + x(1), x(0) * x(1);
    return z;
  }
};

// fun_3_2: R^3 --> R^2 | (x, y, z) --> [(x * y), (y + 2 * z)]
struct fun_3_2 {
  template <typename T>
  inline Matrix<T, Dynamic, 1> operator()(
      const Matrix<T, Dynamic, 1>& x) const {
    Matrix<T, Dynamic, 1> z(2);
    z << x(0) * x(1), x(1) + 2.0 * x(2);
    return z;
  }
};

TEST(FwdFunctor, jacobianMoreOutputsThanInputs) {
  fun_2_3 f;
  Matrix<double, Dynamic, 1> x(2);
  x << 1.5, 2.0;

  Matrix<double, Dynamic, 1> fx;
  Matrix<double, Dynamic, Dynamic> J;
  stan::math::jacobian<double>(f, x, fx, J);

  EXPECT_EQ(3, fx.size());
  EXPECT_FLOAT_EQ(x(0), fx(0));
  EXPECT_FLOAT_EQ(x(0) + x(1), fx(1));
  EXPECT_FLOAT_EQ(x(0) * x(1), fx(2));

  EXPECT_EQ(3, J.rows());
  EXPECT_EQ(2, J.cols());
  EXPECT_FLOAT_EQ(1, J(0, 0));
  EXPECT_FLOAT_EQ(0, J(0, 1));
  EXPECT_FLOAT_EQ(1, J(1, 0));
  EXPECT_FLOAT_EQ(1, J(1, 1));
  EXPECT_FLOAT_EQ(x(1), J(2, 0));
  EXPECT_FLOAT_EQ(x(0), J(2, 1));
}

TEST(FwdFunctor, jacobianFewerOutputsThanInputs) {
  fun_3_2 f;
  Matrix<double, Dynamic, 1> x(3);
  x << 1.5, 2.0, -3.0;

  Matrix<double, Dynamic, 1> fx;
  Matrix<double, Dynamic, Dynamic> J;
  stan::math::jacobian<double>(f, x, fx, J);

  EXPECT_EQ(2, fx.size());
  EXPECT_FLOAT_EQ(x(0) * x(1), fx(0));
  EXPECT_FLOAT_EQ(x(1) + 2.0 * x(2), fx(1));

  EXPECT_EQ(2, J.rows());
  EXPECT_EQ(3, J.cols());
  EXPECT_FLOAT_EQ(x(1), J(0, 0));
  EXPECT_FLOAT_EQ(x(0), J(0, 1));
  EXPECT_FLOAT_EQ(0, J(0, 2));
  EXPECT_FLOAT_EQ(0, J(1, 0));
  EXPECT_FLOAT_EQ(1, J(1, 1));
  EXPECT_FLOAT_EQ(2, J(1, 2));
}
