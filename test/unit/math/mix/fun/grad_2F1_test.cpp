#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(ProbInternalMath, grad2F1_fnegative_z) {
  using stan::math::fvar;

  fvar<double> a = 3.70975;
  fvar<double> b = 1;
  fvar<double> c = 2.70975;
  fvar<double> z = -0.2;

  fvar<double> gradA;
  fvar<double> gradB;
  fvar<double> gradC;
  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);
  EXPECT_NEAR(-0.0488658806159776, gradA.val_, 1e-9);
  EXPECT_NEAR(-0.193844936204681, gradB.val_, 1e-9);
  EXPECT_NEAR(0.0677809985598383, gradC.val_, 1e-9);
}

TEST(ProbInternalMath, grad2F1_fd1) {
  using stan::math::fvar;

  fvar<double> a = 2;
  a.d_ = 1;
  fvar<double> b = 1;
  fvar<double> c = 2;
  fvar<double> z = 0.4;
  fvar<double> gradA;
  fvar<double> gradB;
  fvar<double> gradC;
  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);

  EXPECT_NEAR(0.46177343153972, gradA.val_, 1e-8);
  EXPECT_NEAR(0.16371487651638, gradA.d_, 1e-8);
  EXPECT_NEAR(0.85137603960998, gradB.val_, 1e-8);
  EXPECT_NEAR(-0.4617734352303, gradC.val_, 1e-8);
}
TEST(ProbInternalMath, grad2F1_fd2) {
  using stan::math::fvar;

  fvar<double> a = 2;
  fvar<double> b = 1;
  fvar<double> c = 2;
  b.d_ = 1;
  fvar<double> z = 0.4;
  fvar<double> gradA;
  fvar<double> gradB;
  fvar<double> gradC;
  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);

  EXPECT_NEAR(0.461773431539720, gradA.val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, gradB.val_, 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_, 1e-8);
  EXPECT_NEAR(0.434904696493189, gradB.d_, 1e-8);
}

TEST(ProbInternalMath, grad2F1_fd3) {
  using stan::math::fvar;

  fvar<double> a = 2;
  fvar<double> b = 1;
  fvar<double> c = 2;
  c.d_ = 1;
  fvar<double> z = 0.4;
  fvar<double> gradA;
  fvar<double> gradB;
  fvar<double> gradC;
  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);

  EXPECT_NEAR(0.461773431539720, gradA.val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, gradB.val_, 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_, 1e-8);
  EXPECT_NEAR(0.574406330443730, gradC.d_, 1e-8);
}
TEST(ProbInternalMath, grad2F1_ffd1) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 2;
  a.d_ = 1;
  fvar<fvar<double> > b = 1;
  fvar<fvar<double> > c = 2;
  fvar<fvar<double> > z = 0.4;
  fvar<fvar<double> > gradA;
  fvar<fvar<double> > gradB;
  fvar<fvar<double> > gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val_, 1e-8);
  EXPECT_NEAR(0.163714876516383, gradA.d_.val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, gradB.val_.val_, 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_.val_, 1e-8);
}
TEST(ProbInternalMath, grad2F1_ffd2) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 2.0;
  fvar<fvar<double> > b = 1.0;
  fvar<fvar<double> > c = 2.0;
  b.d_ = 1.0;
  fvar<fvar<double> > z = 0.4;
  fvar<fvar<double> > gradA;
  fvar<fvar<double> > gradB;
  fvar<fvar<double> > gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, gradB.val_.val_, 1e-8);
  EXPECT_NEAR(0.434904696493189, gradB.d_.val_, 1e-8);
  EXPECT_NEAR(-0.461773435230326, gradC.val_.val_, 1e-8);
}

TEST(ProbInternalMath, grad2F1_ffd3) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 2.0;
  fvar<fvar<double> > b = 1.0;
  fvar<fvar<double> > c = 2.0;
  c.d_ = 1.0;
  fvar<fvar<double> > z = 0.4;
  fvar<fvar<double> > gradA;
  fvar<fvar<double> > gradB;
  fvar<fvar<double> > gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, gradB.val_.val_, 1e-8);
  EXPECT_NEAR(-0.461773435230326, gradC.val_.val_, 1e-8);
  EXPECT_NEAR(0.574406330443730, gradC.d_.val_, 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv1) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  a.d_ = 1.0;
  fvar<var> b = 1.0;
  fvar<var> c = 2.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradB;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val(), 1e-8);
  EXPECT_NEAR(0.163714876516383, gradA.d_.val(), 1e-8);
  EXPECT_NEAR(0.851376039609984, gradB.val_.val(), 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_.val(), 1e-8);
}
TEST(ProbInternalMath, grad2F1_fv2) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  fvar<var> b = 1.0;
  b.d_ = 1.0;
  fvar<var> c = 2.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradB;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val(), 1e-8);
  EXPECT_NEAR(0.851376039609984, gradB.val_.val(), 1e-8);
  EXPECT_NEAR(0.434904696493189, gradB.d_.val(), 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_.val(), 1e-8);
}
TEST(ProbInternalMath, grad2F1_fv3) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  fvar<var> b = 1.0;
  fvar<var> c = 2.0;
  c.d_ = 1.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradB;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val(), 1e-8);
  EXPECT_NEAR(0.851376039609984, gradB.val_.val(), 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_.val(), 1e-8);
  EXPECT_NEAR(0.574406330443730, gradC.d_.val(), 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv_1stderiv1) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  a.d_ = 1.0;
  fvar<var> b = 1.0;
  fvar<var> c = 2.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradB;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);

  std::vector<stan::math::var> y1{a.val_};
  std::vector<double> grad1;
  gradA.val_.grad(y1, grad1);
  EXPECT_NEAR(0.163714876516383, grad1[0], 1e-8);
}
TEST(ProbInternalMath, grad2F1_fv_1stderiv2) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  fvar<var> b = 1.0;
  b.d_ = 1.0;
  fvar<var> c = 2.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradB;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);

  std::vector<stan::math::var> y1{b.val_};
  std::vector<double> grad1;
  gradB.val_.grad(y1, grad1);
  EXPECT_NEAR(0.434904696493189, grad1[0], 1e-8);
}
TEST(ProbInternalMath, grad2F1_fv_1stderiv3) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  fvar<var> b = 1.0;
  fvar<var> c = 2.0;
  c.d_ = 1.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradB;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);

  std::vector<stan::math::var> y1{c.val_};
  std::vector<double> grad1;
  gradC.val_.grad(y1, grad1);
  EXPECT_NEAR(0.574406330443730, grad1[0], 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv_2ndderiv1) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  a.d_ = 1.0;
  fvar<var> b = 1.0;
  fvar<var> c = 2.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradB;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);

  std::vector<stan::math::var> y1{a.val_};
  std::vector<double> grad1;
  gradA.d_.grad(y1, grad1);
  EXPECT_NEAR(0.064256527613079, grad1[0], 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv_2ndderiv2) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  fvar<var> b = 1.0;
  b.d_ = 1.0;
  fvar<var> c = 2.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradB;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);

  std::vector<stan::math::var> y1{b.val_};
  std::vector<double> grad1;
  gradB.d_.grad(y1, grad1);
  EXPECT_NEAR(0.222160462864892, grad1[0], 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv_2ndderiv3) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  fvar<var> b = 1.0;
  fvar<var> c = 2.0;
  c.d_ = 1.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradB;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradB, gradC, a, b, c, z);

  std::vector<stan::math::var> y1{c.val_};
  std::vector<double> grad1;
  gradC.d_.grad(y1, grad1);
  EXPECT_NEAR(-1.00024553725447, grad1[0], 1e-8);
}
