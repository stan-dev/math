#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(ProbInternalMath, grad2F1_zero_z) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<double> a1 = 3.70975;
  fvar<double> a2 = 1;
  fvar<double> b1 = 2.70975;
  fvar<double> z = 0;

  a1.d_ = 1;
  a2.d_ = 1;
  b1.d_ = 1;

  auto grad_tuple = stan::math::grad_2F1<true>(a1, a2, b1, z);
  EXPECT_FLOAT_EQ(0, std::get<0>(grad_tuple).val_);
  EXPECT_FLOAT_EQ(0, std::get<1>(grad_tuple).val_);
  EXPECT_FLOAT_EQ(0, std::get<2>(grad_tuple).val_);

  EXPECT_FLOAT_EQ(0, std::get<0>(grad_tuple).d_);
  EXPECT_FLOAT_EQ(0, std::get<1>(grad_tuple).d_);
  EXPECT_FLOAT_EQ(0, std::get<2>(grad_tuple).d_);
}

TEST(ProbInternalMath, grad2F1_fnegative_z) {
  using stan::math::fvar;

  fvar<double> a = 3.70975;
  fvar<double> b = 1;
  fvar<double> c = 2.70975;
  fvar<double> z = -0.2;

  a.d_ = 1;
  b.d_ = 1;
  c.d_ = 1;

  auto grad_tuple = stan::math::grad_2F1(a, b, c, z);
  EXPECT_NEAR(-0.0488658806159776, std::get<0>(grad_tuple), 1e-9);
  EXPECT_NEAR(-0.193844936204681, std::get<1>(grad_tuple), 1e-9);
  EXPECT_NEAR(0.0677809985598383, std::get<2>(grad_tuple), 1e-9);
}

TEST(ProbInternalMath, grad2F1_fd1) {
  using stan::math::fvar;

  fvar<double> a = 2;
  a.d_ = 1;
  fvar<double> b = 1;
  fvar<double> c = 2;
  fvar<double> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);

  EXPECT_NEAR(0.46177343153972, std::get<0>(grad_tuple).val_, 1e-8);
  EXPECT_NEAR(0.16371487651638, std::get<0>(grad_tuple).d_, 1e-8);
  EXPECT_NEAR(0.85137603960998, std::get<1>(grad_tuple).val_, 1e-8);
  EXPECT_NEAR(-0.4617734352303, std::get<2>(grad_tuple).val_, 1e-8);
}
TEST(ProbInternalMath, grad2F1_fd2) {
  using stan::math::fvar;

  fvar<double> a = 2;
  fvar<double> b = 1;
  fvar<double> c = 2;
  b.d_ = 1;
  fvar<double> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);

  EXPECT_NEAR(0.461773431539720, std::get<0>(grad_tuple).val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, std::get<1>(grad_tuple).val_, 1e-8);
  EXPECT_NEAR(-0.46177343523032, std::get<2>(grad_tuple).val_, 1e-8);
  EXPECT_NEAR(0.434904696493189, std::get<1>(grad_tuple).d_, 1e-8);
}

TEST(ProbInternalMath, grad2F1_fd3) {
  using stan::math::fvar;

  fvar<double> a = 2;
  fvar<double> b = 1;
  fvar<double> c = 2;
  c.d_ = 1;
  fvar<double> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);

  EXPECT_NEAR(0.461773431539720, std::get<0>(grad_tuple).val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, std::get<1>(grad_tuple).val_, 1e-8);
  EXPECT_NEAR(-0.46177343523032, std::get<2>(grad_tuple).val_, 1e-8);
  EXPECT_NEAR(0.574406330443730, std::get<2>(grad_tuple).d_, 1e-8);
}
TEST(ProbInternalMath, grad2F1_ffd1) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 2;
  a.d_ = 1;
  fvar<fvar<double> > b = 1;
  fvar<fvar<double> > c = 2;
  fvar<fvar<double> > z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);
  EXPECT_NEAR(0.461773431539720, std::get<0>(grad_tuple).val_.val_, 1e-8);
  EXPECT_NEAR(0.163714876516383, std::get<0>(grad_tuple).d_.val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, std::get<1>(grad_tuple).val_.val_, 1e-8);
  EXPECT_NEAR(-0.46177343523032, std::get<2>(grad_tuple).val_.val_, 1e-8);
}
TEST(ProbInternalMath, grad2F1_ffd2) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 2;
  fvar<fvar<double> > b = 1;
  fvar<fvar<double> > c = 2;
  b.d_ = 1;
  fvar<fvar<double> > z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);
  EXPECT_NEAR(0.461773431539720, std::get<0>(grad_tuple).val_.val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, std::get<1>(grad_tuple).val_.val_, 1e-8);
  EXPECT_NEAR(0.434904696493189, std::get<1>(grad_tuple).d_.val_, 1e-8);
  EXPECT_NEAR(-0.461773435230326, std::get<2>(grad_tuple).val_.val_, 1e-8);
}

TEST(ProbInternalMath, grad2F1_ffd3) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 2;
  fvar<fvar<double> > b = 1;
  fvar<fvar<double> > c = 2;
  c.d_ = 1;
  fvar<fvar<double> > z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);
  EXPECT_NEAR(0.461773431539720, std::get<0>(grad_tuple).val_.val_, 1e-8);
  EXPECT_NEAR(0.851376039609984, std::get<1>(grad_tuple).val_.val_, 1e-8);
  EXPECT_NEAR(-0.461773435230326, std::get<2>(grad_tuple).val_.val_, 1e-8);
  EXPECT_NEAR(0.574406330443730, std::get<2>(grad_tuple).d_.val_, 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv1) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2;
  a.d_ = 1;
  fvar<var> b = 1;
  fvar<var> c = 2;
  fvar<var> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);
  EXPECT_NEAR(0.461773431539720, std::get<0>(grad_tuple).val_.val(), 1e-8);
  EXPECT_NEAR(0.163714876516383, std::get<0>(grad_tuple).d_.val(), 1e-8);
  EXPECT_NEAR(0.851376039609984, std::get<1>(grad_tuple).val_.val(), 1e-8);
  EXPECT_NEAR(-0.46177343523032, std::get<2>(grad_tuple).val_.val(), 1e-8);
}
TEST(ProbInternalMath, grad2F1_fv2) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2;
  fvar<var> b = 1;
  b.d_ = 1;
  fvar<var> c = 2;
  fvar<var> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);
  EXPECT_NEAR(0.461773431539720, std::get<0>(grad_tuple).val_.val(), 1e-8);
  EXPECT_NEAR(0.851376039609984, std::get<1>(grad_tuple).val_.val(), 1e-8);
  EXPECT_NEAR(0.434904696493189, std::get<1>(grad_tuple).d_.val(), 1e-8);
  EXPECT_NEAR(-0.46177343523032, std::get<2>(grad_tuple).val_.val(), 1e-8);
}
TEST(ProbInternalMath, grad2F1_fv3) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2;
  fvar<var> b = 1;
  fvar<var> c = 2;
  c.d_ = 1;
  fvar<var> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);
  EXPECT_NEAR(0.461773431539720, std::get<0>(grad_tuple).val_.val(), 1e-8);
  EXPECT_NEAR(0.851376039609984, std::get<1>(grad_tuple).val_.val(), 1e-8);
  EXPECT_NEAR(-0.46177343523032, std::get<2>(grad_tuple).val_.val(), 1e-8);
  EXPECT_NEAR(0.574406330443730, std::get<2>(grad_tuple).d_.val(), 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv_1stderiv1) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2;
  a.d_ = 1;
  fvar<var> b = 1;
  fvar<var> c = 2;
  fvar<var> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);

  std::vector<stan::math::var> y1{a.val_};
  std::vector<double> grad1;
  std::get<0>(grad_tuple).val_.grad(y1, grad1);
  EXPECT_NEAR(0.163714876516383, grad1[0], 1e-8);
}
TEST(ProbInternalMath, grad2F1_fv_1stderiv2) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2;
  fvar<var> b = 1;
  b.d_ = 1;
  fvar<var> c = 2;
  fvar<var> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);

  std::vector<stan::math::var> y1{b.val_};
  std::vector<double> grad1;
  std::get<1>(grad_tuple).val_.grad(y1, grad1);
  EXPECT_NEAR(0.434904696493189, grad1[0], 1e-8);
}
TEST(ProbInternalMath, grad2F1_fv_1stderiv3) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2;
  fvar<var> b = 1;
  fvar<var> c = 2;
  c.d_ = 1;
  fvar<var> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);

  std::vector<stan::math::var> y1{c.val_};
  std::vector<double> grad1;
  std::get<2>(grad_tuple).val_.grad(y1, grad1);
  EXPECT_NEAR(0.574406330443730, grad1[0], 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv_2ndderiv1) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2;
  a.d_ = 1;
  fvar<var> b = 1;
  fvar<var> c = 2;
  fvar<var> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);

  std::vector<stan::math::var> y1{a.val_};
  std::vector<double> grad1;
  std::get<0>(grad_tuple).d_.grad(y1, grad1);
  EXPECT_NEAR(0.064256527613079, grad1[0], 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv_2ndderiv2) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2;
  fvar<var> b = 1;
  b.d_ = 1;
  fvar<var> c = 2;
  fvar<var> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);

  std::vector<stan::math::var> y1{b.val_};
  std::vector<double> grad1;
  std::get<1>(grad_tuple).d_.grad(y1, grad1);
  EXPECT_NEAR(0.222160462864892, grad1[0], 1e-8);
}

TEST(ProbInternalMath, grad2F1_fv_2ndderiv3) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 2;
  fvar<var> b = 1;
  fvar<var> c = 2;
  c.d_ = 1;
  fvar<var> z = 0.4;

  auto grad_tuple = stan::math::grad_2F1<true>(a, b, c, z);

  std::vector<stan::math::var> y1{c.val_};
  std::vector<double> grad1;
  std::get<2>(grad_tuple).d_.grad(y1, grad1);
  EXPECT_NEAR(-1.00024553725447, grad1[0], 1e-8);
}
