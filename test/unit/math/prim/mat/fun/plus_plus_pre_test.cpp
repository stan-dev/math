#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>


template <typename C>
void test_val_ref(C& x) {
  using stan::math::plus_plus_pre;
  // can't use << because it's row-major and need col-major for matrix
  x(0) = 1;
  x(1) = 2.5;
  x(2) = -1.4;
  x(3) = -3;
  C y = plus_plus_pre(x);

  // test self update
  EXPECT_FLOAT_EQ(2, x(0));
  EXPECT_FLOAT_EQ(3.5, x(1));
  EXPECT_FLOAT_EQ(-0.4, x(2));

  // test return
  EXPECT_FLOAT_EQ(2, y(0));
  EXPECT_FLOAT_EQ(3.5, y(1));
  EXPECT_FLOAT_EQ(-0.4, y(2));
  
  // test reference
  EXPECT_EQ(&x, &plus_plus_pre(x));
}

TEST(MathFunctions, plusPlusPreVec) {
  Eigen::VectorXd v(4);
  test_val_ref(v);
}
TEST(MathFunctions, plusPlusPreRowVec) {
  Eigen::RowVectorXd rv(4);
  test_val_ref(rv);
}
TEST(MathFunctions, plusPlusPreMat) {
  Eigen::MatrixXd m(2, 2);
  test_val_ref(m);
}
