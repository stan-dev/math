#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <vector>

TEST(MathPrimScalFun, grad2F1_negative_z) {
  double a1 = 3.70975;
  double a2 = 1;
  double b1 = 2.70975;
  double z = -0.2;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_NEAR(-0.0488658806159776, std::get<0>(grad_tuple), 1e-9);
  EXPECT_NEAR(-0.193844936204681, std::get<1>(grad_tuple), 1e-9);
  EXPECT_NEAR(0.0677809985598383, std::get<2>(grad_tuple), 1e-9);
}

TEST(MathPrimScalFun, grad2F1_zero_z) {
  double a1 = 3.70975;
  double a2 = 1;
  double b1 = 2.70975;
  double z = 0;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_FLOAT_EQ(0, std::get<0>(grad_tuple));
  EXPECT_FLOAT_EQ(0, std::get<1>(grad_tuple));
  EXPECT_FLOAT_EQ(0, std::get<2>(grad_tuple));
}

TEST(MathPrimScalFun, grad2F1_1) {
  double a1 = 1;
  double a2 = 1;
  double b1 = 1;
  double z = 0.6;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_NEAR(2.290726829685388, std::get<0>(grad_tuple), 1e-9);
  EXPECT_NEAR(2.290726829685388, std::get<1>(grad_tuple), 1e-9);
  EXPECT_NEAR(-2.290726829685388, std::get<2>(grad_tuple), 1e-9);
}

TEST(MathPrimScalFun, grad2F1_2) {
  double a1 = 1;
  double a2 = 31;
  double b1 = 41;
  double z = 1;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_NEAR(6.825270649241036, std::get<0>(grad_tuple), 1e-8);
  EXPECT_NEAR(0.4938271604938271, std::get<1>(grad_tuple), 1e-8);
  EXPECT_NEAR(-0.382716049382716, std::get<2>(grad_tuple), 1e-8);
}

TEST(MathPrimScalFun, grad2F1_3) {
  double a1 = 1;
  double a2 = -2.1;
  double b1 = 41;
  double z = 1;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_NEAR(-0.04921317604093563, std::get<0>(grad_tuple), 1e-8);
  EXPECT_NEAR(0.02256814168279349, std::get<1>(grad_tuple), 1e-8);
  EXPECT_NEAR(0.00118482743834665, std::get<2>(grad_tuple), 1e-8);
}

TEST(MathPrimScalFun, grad2F1_4) {
  double a1 = 1;
  double a2 = 12;
  double b1 = 10;
  double z = 1;

  EXPECT_THROW(auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z),
               std::domain_error);
}

TEST(MathPrimScalFun, grad2F1_5) {
  double a1 = 1;
  double a2 = 12;
  double b1 = 20;
  double z = 1.2;

  EXPECT_THROW(auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z),
               std::domain_error);
}

TEST(MathPrimScalFun, grad2F1_6) {
  double a1 = 1;
  double a2 = -0.5;
  double b1 = 10.6;
  double z = 0.3;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_NEAR(-0.01443822031245647, std::get<0>(grad_tuple), 1e-8);
  EXPECT_NEAR(0.02829710651967078, std::get<1>(grad_tuple), 1e-8);
  EXPECT_NEAR(0.00136986255602642, std::get<2>(grad_tuple), 1e-8);
}

TEST(MathPrimScalFun, grad2F1_7) {
  double a1 = 1;
  double a2 = -0.5;
  double b1 = 10;
  double z = 0.3;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_NEAR(-0.0153218866216130, std::get<0>(grad_tuple), 1e-8);
  EXPECT_NEAR(0.02999436412836072, std::get<1>(grad_tuple), 1e-8);
  EXPECT_NEAR(0.0015413242328729, std::get<2>(grad_tuple), 1e-8);
}

TEST(MathPrimScalFun, grad2F1_8) {
  double a1 = -.5;
  double a2 = -4.5;
  double b1 = 11;
  double z = 0.3;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_NEAR(-0.1227022810085707, std::get<0>(grad_tuple), 1e-8);
  EXPECT_NEAR(-0.01298849638043795, std::get<1>(grad_tuple), 1e-8);
  EXPECT_NEAR(-0.0053540982315572, std::get<2>(grad_tuple), 1e-8);
}

TEST(MathPrimScalFun, grad2F1_9) {
  double a1 = -.5;
  double a2 = -4.5;
  double b1 = -3.2;
  double z = 0.9;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_NEAR(0.85880025358111, std::get<0>(grad_tuple), 1e-8);
  EXPECT_NEAR(0.4677704416159314, std::get<1>(grad_tuple), 1e-8);
  EXPECT_NEAR(-4.19010422485256, std::get<2>(grad_tuple), 1e-8);
}

TEST(MathPrimScalFun, grad2F1_10) {
  double a1 = 2;
  double a2 = 1;
  double b1 = 2;
  double z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_NEAR(0.4617734323582945, std::get<0>(grad_tuple), 1e-8);
  EXPECT_NEAR(0.851376039609984, std::get<1>(grad_tuple), 1e-8);
  EXPECT_NEAR(-0.4617734323582945, std::get<2>(grad_tuple), 1e-8);
}

TEST(MathPrimScalFun, grad2F1_11) {
  double a1 = 3.70975;
  double a2 = 1;
  double b1 = 2.70975;
  double z = 0.999696;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);

  EXPECT_NEAR(29369830.002773938200417693317785, std::get<0>(grad_tuple),
              1e-1);  // reference: discrete diff in mathematica
  EXPECT_NEAR(36347869.41885337, std::get<1>(grad_tuple), 1e-1);
  EXPECT_NEAR(-30843032.10697079073015067426929807, std::get<2>(grad_tuple),
              1e-1);  // reference: discrete diff in mathematica
}
