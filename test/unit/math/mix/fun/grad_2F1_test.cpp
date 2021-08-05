#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(ProbInternalMath, grad2F1_fd1) {
  using stan::math::fvar;

  fvar<double> a = 2.0;
  a.d_ = 1.0;
  fvar<double> b = 1.0;
  fvar<double> c = 2.0;
  fvar<double> z = 0.4;
  fvar<double> gradA;
  fvar<double> gradC;
  stan::math::grad_2F1(gradA, gradC, a, b, c, z);

  EXPECT_NEAR(0.461773431539720, gradA.val_, 1e-8);
  EXPECT_NEAR(0.163714876516383, gradA.d_, 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_, 1e-8);
}
TEST(ProbInternalMath, grad2F1_fd2) {
  using stan::math::fvar;

  fvar<double> a = 2.0;
  fvar<double> b = 1.0;
  fvar<double> c = 2.0;
  c.d_ = 1.0;
  fvar<double> z = 0.4;
  fvar<double> gradA;
  fvar<double> gradC;
  stan::math::grad_2F1(gradA, gradC, a, b, c, z);

  EXPECT_NEAR(0.461773431539720, gradA.val_, 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_, 1e-8);
  EXPECT_NEAR(0.574406330443730, gradC.d_, 1e-8);
}
TEST(ProbInternalMath, grad2F1_ffd1) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 2.0;
  a.d_ = 1.0;
  fvar<fvar<double> > b = 1.0;
  fvar<fvar<double> > c = 2.0;
  fvar<fvar<double> > z = 0.4;
  fvar<fvar<double> > gradA;
  fvar<fvar<double> > gradC;

  stan::math::grad_2F1(gradA, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val_, 1e-8);
  EXPECT_NEAR(0.163714876516383, gradA.d_.val_, 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_.val_, 1e-8);
}
TEST(ProbInternalMath, grad2F1_ffd2) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 2.0;
  fvar<fvar<double> > b = 1.0;
  fvar<fvar<double> > c = 2.0;
  c.d_ = 1.0;
  fvar<fvar<double> > z = 0.4;
  fvar<fvar<double> > gradA;
  fvar<fvar<double> > gradC;

  stan::math::grad_2F1(gradA, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val_, 1e-8);
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
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val(), 1e-8);
  EXPECT_NEAR(0.163714876516383, gradA.d_.val(), 1e-8);
  EXPECT_NEAR(-0.46177343523032, gradC.val_.val(), 1e-8);
}
TEST(ProbInternalMath, grad2F1_fv2) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2.0;
  fvar<var> b = 1.0;
  fvar<var> c = 2.0;
  c.d_ = 1.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradC, a, b, c, z);
  EXPECT_NEAR(0.461773431539720, gradA.val_.val(), 1e-8);
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
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradC, a, b, c, z);

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
  fvar<var> c = 2.0;
  c.d_ = 1.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradC, a, b, c, z);

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
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradC, a, b, c, z);

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
  fvar<var> c = 2.0;
  c.d_ = 1.0;
  fvar<var> z = 0.4;
  fvar<var> gradA;
  fvar<var> gradC;

  stan::math::grad_2F1(gradA, gradC, a, b, c, z);

  std::vector<stan::math::var> y1{c.val_};
  std::vector<double> grad1;
  gradC.d_.grad(y1, grad1);
  EXPECT_NEAR(-1.00024553725447, grad1[0], 1e-8);
}
