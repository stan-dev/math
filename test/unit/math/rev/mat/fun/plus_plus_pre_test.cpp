#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

template <typename C>
void expect_vals_derivs(C& x) {
  using stan::math::plus_plus_pre;
  // can't use << because it's row-major and need col-major for matrix
  x(0) = 1;
  x(1) = 2.5;
  x(2) = -1.4;
  x(3) = -3;
  C y = plus_plus_pre(x);

  // test self update
  EXPECT_FLOAT_EQ(2, x(0).val());
  EXPECT_FLOAT_EQ(3.5, x(1).val());
  EXPECT_FLOAT_EQ(-0.4, x(2).val());

  // test return
  EXPECT_FLOAT_EQ(2, y(0).val());
  EXPECT_FLOAT_EQ(3.5, y(1).val());
  EXPECT_FLOAT_EQ(-0.4, y(2).val());
  
  // test derivative
  AVEC xs = createAVEC(x(0), x(1), x(2));
  VEC g;
  // just choose an output element to test
  y(1).grad(xs, g);
  EXPECT_EQ(3, g.size());
  EXPECT_FLOAT_EQ(0, g[0]);
  EXPECT_FLOAT_EQ(1, g[1]);
  EXPECT_FLOAT_EQ(0, g[2]);

  // test reference
  EXPECT_EQ(&x, &plus_plus_pre(x));
}

TEST(MathFunctions, plusPlusPreVec) {
  using stan::math::var;
  Eigen::Matrix<var, -1, 1> v(4);
  expect_vals_derivs(v);
}
TEST(MathFunctions, plusPlusPreRowVec) {
  using stan::math::var;
  Eigen::Matrix<var, 1, -1> rv(4);
  expect_vals_derivs(rv);
}
TEST(MathFunctions, plusPlusPreMat) {
  using stan::math::var;
  Eigen::Matrix<var, -1, -1> m(2, 2);
  expect_vals_derivs(m);
}
