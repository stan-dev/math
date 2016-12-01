#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>


template <typename C>
void test_val_ref(C& x) {
  using stan::math::plus_plus_pre;
  // can't use << because it's row-major and need col-major for matrix
  x(0) = 1;
  x(1) = 2.5;
  x(2) = -1.4;
  x(3) = -3;
  for (int i = 0; i < x.size(); ++i)
    x(i).d_ = 8.5 + i;
  C y = plus_plus_pre(x);

  // test self update
  EXPECT_FLOAT_EQ(2, x(0).val_);
  EXPECT_FLOAT_EQ(3.5, x(1).val_);
  EXPECT_FLOAT_EQ(-0.4, x(2).val_);

  // test return
  EXPECT_FLOAT_EQ(2, y(0).val_);
  EXPECT_FLOAT_EQ(3.5, y(1).val_);
  EXPECT_FLOAT_EQ(-0.4, y(2).val_);

  // test derivative
  EXPECT_FLOAT_EQ(8.5, y(0).d_);
  EXPECT_FLOAT_EQ(9.5, y(1).d_);
  EXPECT_FLOAT_EQ(10.5, y(2).d_);
  
  // test reference
  EXPECT_EQ(&x, &plus_plus_pre(x));
}

TEST(MathFunctions, plusPlusPreVec) {
  using stan::math::fvar;
  Eigen::Matrix<fvar<double>, -1, 1> v(4);
  test_val_ref(v);
}
TEST(MathFunctions, plusPlusPreRowVec) {
  using stan::math::fvar;
  Eigen::Matrix<fvar<double>, 1, -1> rv(4);
  test_val_ref(rv);
}
TEST(MathFunctions, plusPlusPreMat) {
  using stan::math::fvar;
  Eigen::Matrix<fvar<double>, -1, -1> m(2, 2);
  test_val_ref(m);
}

template <typename C>
void test_val_ref2(C& x) {
  using stan::math::plus_plus_pre;
  // can't use << because it's row-major and need col-major for matrix
  x(0) = 1;
  x(1) = 2.5;
  x(2) = -1.4;
  x(3) = -3;
  for (int i = 0; i < x.size(); ++i)
    x(i).d_ = 8.5 + i;
  C y = plus_plus_pre(x);

  // test self update
  EXPECT_FLOAT_EQ(2, x(0).val_.val_);
  EXPECT_FLOAT_EQ(3.5, x(1).val_.val_);
  EXPECT_FLOAT_EQ(-0.4, x(2).val_.val_);

  // test return
  EXPECT_FLOAT_EQ(2, y(0).val_.val_);
  EXPECT_FLOAT_EQ(3.5, y(1).val_.val_);
  EXPECT_FLOAT_EQ(-0.4, y(2).val_.val_);

  // test derivative
  EXPECT_FLOAT_EQ(8.5, y(0).d_.val_);
  EXPECT_FLOAT_EQ(9.5, y(1).d_.val_);
  EXPECT_FLOAT_EQ(10.5, y(2).d_.val_);
  
  // test reference
  EXPECT_EQ(&x, &plus_plus_pre(x));
}

TEST(MathFunctions, plusPlusPreVec2) {
  using stan::math::fvar;
  Eigen::Matrix<fvar<fvar<double> >, -1, 1> v(4);
  test_val_ref2(v);
}
TEST(MathFunctions, plusPlusPreRowVec2) {
  using stan::math::fvar;
  Eigen::Matrix<fvar<fvar<double> >, 1, -1> rv(4);
  test_val_ref2(rv);
}
TEST(MathFunctions, plusPlusPreMat2) {
  using stan::math::fvar;
  Eigen::Matrix<fvar<fvar<double> >, -1, -1> m(2, 2);
  test_val_ref2(m);
}
