#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <vector>

TEST(MathPrimScalFun, grad2F1_1) {
  double a1 = 1.0;
  double a2 = 1.0;
  double b1 = 1.0;
  double z = 0.6;

  double grad_a1;
  double grad_b1;
  stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z);
  EXPECT_NEAR(2.290726829685388, grad_a1, 1e-8);
  EXPECT_NEAR(-2.290726829685388, grad_b1, 1e-8);
}

TEST(MathPrimScalFun, grad2F1_2) {
  double a1 = 1.0;
  double a2 = 31.0;
  double b1 = 41.0;
  double z = 1.0;

  double grad_a1;
  double grad_b1;
  stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z, 1e-11);
  EXPECT_NEAR(6.825270649241036, grad_a1, 1e-8);
  EXPECT_NEAR(-0.382716049382716, grad_b1, 1e-8);
}

TEST(MathPrimScalFun, grad2F1_3) {
  double a1 = 1.0;
  double a2 = -2.1;
  double b1 = 41.0;
  double z = 1.0;

  double grad_a1;
  double grad_b1;
  stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z);
  EXPECT_NEAR(-0.04921317604093563, grad_a1, 1e-8);
  EXPECT_NEAR(0.00118482743834665, grad_b1, 1e-8);
}

TEST(MathPrimScalFun, grad2F1_4) {
  double a1 = 1.0;
  double a2 = 12.0;
  double b1 = 10.0;
  double z = 1.0;

  double grad_a1;
  double grad_b1;
  EXPECT_THROW(stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z),
               std::domain_error);
}

TEST(MathPrimScalFun, grad2F1_5) {
  double a1 = 1.0;
  double a2 = 12.0;
  double b1 = 20.0;
  double z = 1.2;

  double grad_a1;
  double grad_b1;
  EXPECT_THROW(stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z),
               std::domain_error);
}

TEST(MathPrimScalFun, grad2F1_6) {
  double a1 = 1.0;
  double a2 = -0.5;
  double b1 = 10.6;
  double z = 0.3;

  double grad_a1;
  double grad_b1;
  stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z);
  EXPECT_NEAR(-0.01443822031245647, grad_a1, 1e-8);
  EXPECT_NEAR(0.00136986255602642, grad_b1, 1e-8);
}

TEST(MathPrimScalFun, grad2F1_7) {
  double a1 = 1.0;
  double a2 = -0.5;
  double b1 = 10.0;
  double z = 0.3;

  double grad_a1;
  double grad_b1;
  stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z);
  EXPECT_NEAR(-0.0153218866216130, grad_a1, 1e-8);
  EXPECT_NEAR(0.0015413242328729, grad_b1, 1e-8);
}

TEST(MathPrimScalFun, grad2F1_8) {
  double a1 = -.5;
  double a2 = -4.5;
  double b1 = 11.0;
  double z = 0.3;

  double grad_a1;
  double grad_b1;
  stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z);
  EXPECT_NEAR(-0.1227022810085707, grad_a1, 1e-8);
  EXPECT_NEAR(-0.0053540982315572, grad_b1, 1e-8);
}

TEST(MathPrimScalFun, grad2F1_9) {
  double a1 = -.5;
  double a2 = -4.5;
  double b1 = -3.2;
  double z = 0.9;

  double grad_a1;
  double grad_b1;
  stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z);
  EXPECT_NEAR(0.85880025358111, grad_a1, 1e-8);
  EXPECT_NEAR(-4.19010422485256, grad_b1, 1e-8);
}

TEST(MathPrimScalFun, grad2F1_10) {
  double a1 = 2;
  double a2 = 1.0;
  double b1 = 2;
  double z = 0.4;

  double grad_a1;
  double grad_b1;
  stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z);
  EXPECT_NEAR(0.4617734323582945, grad_a1, 1e-8);
  EXPECT_NEAR(-0.4617734323582945, grad_b1, 1e-8);
}

TEST(MathPrimScalFun, grad2F1_11) {
  double a1 = 3.70975;
  double a2 = 1.0;
  double b1 = 2.70975;
  double z = 0.999696;

  double grad_a1;
  double grad_b1;
  stan::math::grad_2F1(grad_a1, grad_b1, a1, a2, b1, z);

  EXPECT_NEAR(29369830.002773938200417693317785, grad_a1,
              1e-1);  // reference: discrete diff in mathematica
  EXPECT_NEAR(-30843032.10697079073015067426929807, grad_b1,
              1e-1);  // reference: discrete diff in mathematica
}
