#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <complex>

TEST(PrimScalarSigTests, abs) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::abs(real_1);
  auto result_2 = stan::math::abs(int_1);
  auto result_3 = stan::math::abs(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::abs(real_1));
  EXPECT_FLOAT_EQ(result_2, std::abs(int_1));
  EXPECT_FLOAT_EQ(result_3, std::abs(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, acos) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 1;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::acos(real_1);
  auto result_2 = stan::math::acos(int_1);
  auto result_3 = stan::math::acos(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::acos(real_1));
  EXPECT_FLOAT_EQ(result_2, std::acos(int_1));
  EXPECT_EQ(result_3, std::acos(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, acosh) {
  using namespace std::complex_literals;
  double real_1 = 1.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::acosh(real_1);
  auto result_2 = stan::math::acosh(int_1);
  auto result_3 = stan::math::acosh(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::acosh(real_1));
  EXPECT_FLOAT_EQ(result_2, std::acosh(int_1));
  EXPECT_EQ(result_3, std::acosh(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, arg) {
  using namespace std::complex_literals;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::arg(complex_1);
  EXPECT_EQ(result_1, std::arg(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, asin) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 1;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::asin(real_1);
  auto result_2 = stan::math::asin(int_1);
  auto result_3 = stan::math::asin(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::asin(real_1));
  EXPECT_FLOAT_EQ(result_2, std::asin(int_1));
  EXPECT_EQ(result_3, std::asin(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, asinh) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::asinh(real_1);
  auto result_2 = stan::math::asinh(int_1);
  auto result_3 = stan::math::asinh(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::asinh(real_1));
  EXPECT_FLOAT_EQ(result_2, std::asinh(int_1));
  EXPECT_EQ(result_3, std::asinh(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, atan) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::atan(real_1);
  auto result_2 = stan::math::atan(int_1);
  auto result_3 = stan::math::atan(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::atan(real_1));
  EXPECT_FLOAT_EQ(result_2, std::atan(int_1));
  EXPECT_EQ(result_3, std::atan(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, atanh) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 1;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::atanh(real_1);
  auto result_2 = stan::math::atanh(int_1);
  auto result_3 = stan::math::atanh(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::atanh(real_1));
  EXPECT_FLOAT_EQ(result_2, std::atanh(int_1));
  EXPECT_EQ(result_3, std::atanh(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, cbrt) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::cbrt(real_1);
  auto result_2 = stan::math::cbrt(int_1);
  EXPECT_FLOAT_EQ(result_1, std::cbrt(real_1));
  EXPECT_FLOAT_EQ(result_2, std::cbrt(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, ceil) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::ceil(real_1);
  auto result_2 = stan::math::ceil(int_1);
  EXPECT_FLOAT_EQ(result_1, std::ceil(real_1));
  EXPECT_FLOAT_EQ(result_2, std::ceil(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, conj) {
  using namespace std::complex_literals;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::conj(complex_1);
  EXPECT_EQ(result_1, std::conj(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, cos) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::cos(real_1);
  auto result_2 = stan::math::cos(int_1);
  auto result_3 = stan::math::cos(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::cos(real_1));
  EXPECT_FLOAT_EQ(result_2, std::cos(int_1));
  EXPECT_EQ(result_3, std::cos(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, cosh) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::cosh(real_1);
  auto result_2 = stan::math::cosh(int_1);
  auto result_3 = stan::math::cosh(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::cosh(real_1));
  EXPECT_FLOAT_EQ(result_2, std::cosh(int_1));
  EXPECT_EQ(result_3, std::cosh(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, erf) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::erf(real_1);
  auto result_2 = stan::math::erf(int_1);
  EXPECT_FLOAT_EQ(result_1, std::erf(real_1));
  EXPECT_FLOAT_EQ(result_2, std::erf(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, erfc) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::erfc(real_1);
  auto result_2 = stan::math::erfc(int_1);
  EXPECT_FLOAT_EQ(result_1, std::erfc(real_1));
  EXPECT_FLOAT_EQ(result_2, std::erfc(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, exp) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::exp(real_1);
  auto result_2 = stan::math::exp(int_1);
  auto result_3 = stan::math::exp(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::exp(real_1));
  EXPECT_FLOAT_EQ(result_2, std::exp(int_1));
  EXPECT_EQ(result_3, std::exp(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, exp2) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::exp2(real_1);
  auto result_2 = stan::math::exp2(int_1);
  EXPECT_FLOAT_EQ(result_1, std::exp2(real_1));
  EXPECT_FLOAT_EQ(result_2, std::exp2(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, expm1) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::expm1(real_1);
  auto result_2 = stan::math::expm1(int_1);
  EXPECT_FLOAT_EQ(result_1, std::expm1(real_1));
  EXPECT_FLOAT_EQ(result_2, std::expm1(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, fabs) {
  double real_1 = -0.4;
  int int_1 = 2;
  auto result_1 = stan::math::fabs(real_1);
  auto result_2 = stan::math::fabs(int_1);
  EXPECT_FLOAT_EQ(result_1, std::fabs(real_1));
  EXPECT_FLOAT_EQ(result_2, std::fabs(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, floor) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::floor(real_1);
  auto result_2 = stan::math::floor(int_1);
  EXPECT_FLOAT_EQ(result_1, std::floor(real_1));
  EXPECT_FLOAT_EQ(result_2, std::floor(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, lgamma) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::lgamma(real_1);
  auto result_2 = stan::math::lgamma(int_1);
  EXPECT_FLOAT_EQ(result_1, std::lgamma(real_1));
  EXPECT_FLOAT_EQ(result_2, std::lgamma(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, log) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::log(real_1);
  auto result_2 = stan::math::log(int_1);
  auto result_3 = stan::math::log(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::log(real_1));
  EXPECT_FLOAT_EQ(result_2, std::log(int_1));
  EXPECT_EQ(result_3, std::log(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, log10) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::log10(real_1);
  auto result_2 = stan::math::log10(int_1);
  auto result_3 = stan::math::log10(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::log10(real_1));
  EXPECT_FLOAT_EQ(result_2, std::log10(int_1));
  EXPECT_EQ(result_3, std::log10(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, log1p) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::log1p(real_1);
  auto result_2 = stan::math::log1p(int_1);
  EXPECT_FLOAT_EQ(result_1, std::log1p(real_1));
  EXPECT_FLOAT_EQ(result_2, std::log1p(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, log2) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::log2(real_1);
  auto result_2 = stan::math::log2(int_1);
  EXPECT_FLOAT_EQ(result_1, std::log2(real_1));
  EXPECT_FLOAT_EQ(result_2, std::log2(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, norm) {
  using namespace std::complex_literals;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::norm(complex_1);
  EXPECT_EQ(result_1, std::norm(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, proj) {
  using namespace std::complex_literals;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::proj(complex_1);
  EXPECT_EQ(result_1, std::proj(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, round) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::round(real_1);
  auto result_2 = stan::math::round(int_1);
  EXPECT_FLOAT_EQ(result_1, std::round(real_1));
  EXPECT_FLOAT_EQ(result_2, std::round(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, sin) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::sin(real_1);
  auto result_2 = stan::math::sin(int_1);
  auto result_3 = stan::math::sin(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::sin(real_1));
  EXPECT_FLOAT_EQ(result_2, std::sin(int_1));
  EXPECT_EQ(result_3, std::sin(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, sinh) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::sinh(real_1);
  auto result_2 = stan::math::sinh(int_1);
  auto result_3 = stan::math::sinh(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::sinh(real_1));
  EXPECT_FLOAT_EQ(result_2, std::sinh(int_1));
  EXPECT_EQ(result_3, std::sinh(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, sqrt) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::sqrt(real_1);
  auto result_2 = stan::math::sqrt(int_1);
  auto result_3 = stan::math::sqrt(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::sqrt(real_1));
  EXPECT_FLOAT_EQ(result_2, std::sqrt(int_1));
  EXPECT_EQ(result_3, std::sqrt(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, tan) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::tan(real_1);
  auto result_2 = stan::math::tan(int_1);
  auto result_3 = stan::math::tan(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::tan(real_1));
  EXPECT_FLOAT_EQ(result_2, std::tan(int_1));
  EXPECT_EQ(result_3, std::tan(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, tanh) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  int int_1 = 2;
  std::complex<double> complex_1 = 1. + 2.0i;
  auto result_1 = stan::math::tanh(real_1);
  auto result_2 = stan::math::tanh(int_1);
  auto result_3 = stan::math::tanh(complex_1);
  EXPECT_FLOAT_EQ(result_1, std::tanh(real_1));
  EXPECT_FLOAT_EQ(result_2, std::tanh(int_1));
  EXPECT_EQ(result_3, std::tanh(complex_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, tgamma) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::tgamma(real_1);
  auto result_2 = stan::math::tgamma(int_1);
  EXPECT_FLOAT_EQ(result_1, std::tgamma(real_1));
  EXPECT_FLOAT_EQ(result_2, std::tgamma(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, trunc) {
  double real_1 = 0.4;
  int int_1 = 2;
  auto result_1 = stan::math::trunc(real_1);
  auto result_2 = stan::math::trunc(int_1);
  EXPECT_FLOAT_EQ(result_1, std::trunc(real_1));
  EXPECT_FLOAT_EQ(result_2, std::trunc(int_1));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, atan2) {
  double real_1 = 0.4;
  double real_2 = 2.2;
  int int_1 = 2;
  int int_2 = 4;
  auto result_1 = stan::math::atan2(real_1, real_2);
  auto result_2 = stan::math::atan2(real_1, int_2);
  auto result_3 = stan::math::atan2(int_1, real_2);
  auto result_4 = stan::math::atan2(int_1, int_2);
  EXPECT_FLOAT_EQ(result_1, std::atan2(real_1, real_2));
  EXPECT_FLOAT_EQ(result_2, std::atan2(real_1, int_2));
  EXPECT_FLOAT_EQ(result_3, std::atan2(int_1, real_2));
  EXPECT_FLOAT_EQ(result_4, std::atan2(int_1, int_2));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, copysign) {
  double real_1 = 0.4;
  double real_2 = 2.2;
  int int_1 = 2;
  int int_2 = 4;
  auto result_1 = stan::math::copysign(real_1, real_2);
  auto result_2 = stan::math::copysign(real_1, int_2);
  auto result_3 = stan::math::copysign(int_1, real_2);
  auto result_4 = stan::math::copysign(int_1, int_2);
  EXPECT_FLOAT_EQ(result_1, std::copysign(real_1, real_2));
  EXPECT_FLOAT_EQ(result_2, std::copysign(real_1, int_2));
  EXPECT_FLOAT_EQ(result_3, std::copysign(int_1, real_2));
  EXPECT_FLOAT_EQ(result_4, std::copysign(int_1, int_2));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, fdim) {
  double real_1 = 0.4;
  double real_2 = 2.2;
  int int_1 = 2;
  int int_2 = 4;
  auto result_1 = stan::math::fdim(real_1, real_2);
  auto result_2 = stan::math::fdim(real_1, int_2);
  auto result_3 = stan::math::fdim(int_1, real_2);
  auto result_4 = stan::math::fdim(int_1, int_2);
  EXPECT_FLOAT_EQ(result_1, std::fdim(real_1, real_2));
  EXPECT_FLOAT_EQ(result_2, std::fdim(real_1, int_2));
  EXPECT_FLOAT_EQ(result_3, std::fdim(int_1, real_2));
  EXPECT_FLOAT_EQ(result_4, std::fdim(int_1, int_2));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, fmax) {
  double real_1 = 0.4;
  double real_2 = 2.2;
  int int_1 = 2;
  int int_2 = 4;
  auto result_1 = stan::math::fmax(real_1, real_2);
  auto result_2 = stan::math::fmax(real_1, int_2);
  auto result_3 = stan::math::fmax(int_1, real_2);
  auto result_4 = stan::math::fmax(int_1, int_2);
  EXPECT_FLOAT_EQ(result_1, std::fmax(real_1, real_2));
  EXPECT_FLOAT_EQ(result_2, std::fmax(real_1, int_2));
  EXPECT_FLOAT_EQ(result_3, std::fmax(int_1, real_2));
  EXPECT_FLOAT_EQ(result_4, std::fmax(int_1, int_2));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, fmin) {
  double real_1 = 0.4;
  double real_2 = 2.2;
  int int_1 = 2;
  int int_2 = 4;
  auto result_1 = stan::math::fmin(real_1, real_2);
  auto result_2 = stan::math::fmin(real_1, int_2);
  auto result_3 = stan::math::fmin(int_1, real_2);
  auto result_4 = stan::math::fmin(int_1, int_2);
  EXPECT_FLOAT_EQ(result_1, std::fmin(real_1, real_2));
  EXPECT_FLOAT_EQ(result_2, std::fmin(real_1, int_2));
  EXPECT_FLOAT_EQ(result_3, std::fmin(int_1, real_2));
  EXPECT_FLOAT_EQ(result_4, std::fmin(int_1, int_2));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, fmod) {
  double real_1 = 0.4;
  double real_2 = 2.2;
  int int_1 = 2;
  int int_2 = 4;
  auto result_1 = stan::math::fmod(real_1, real_2);
  auto result_2 = stan::math::fmod(real_1, int_2);
  auto result_3 = stan::math::fmod(int_1, real_2);
  auto result_4 = stan::math::fmod(int_1, int_2);
  EXPECT_FLOAT_EQ(result_1, std::fmod(real_1, real_2));
  EXPECT_FLOAT_EQ(result_2, std::fmod(real_1, int_2));
  EXPECT_FLOAT_EQ(result_3, std::fmod(int_1, real_2));
  EXPECT_FLOAT_EQ(result_4, std::fmod(int_1, int_2));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, hypot) {
  double real_1 = 0.4;
  double real_2 = 2.2;
  int int_1 = 2;
  int int_2 = 4;
  auto result_1 = stan::math::hypot(real_1, real_2);
  auto result_2 = stan::math::hypot(real_1, int_2);
  auto result_3 = stan::math::hypot(int_1, real_2);
  auto result_4 = stan::math::hypot(int_1, int_2);
  EXPECT_FLOAT_EQ(result_1, std::hypot(real_1, real_2));
  EXPECT_FLOAT_EQ(result_2, std::hypot(real_1, int_2));
  EXPECT_FLOAT_EQ(result_3, std::hypot(int_1, real_2));
  EXPECT_FLOAT_EQ(result_4, std::hypot(int_1, int_2));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, ldexp) {
  double real_1 = 0.4;
  double real_2 = 2.2;
  int int_1 = 2;
  int int_2 = 4;
  auto result_1 = stan::math::ldexp(real_1, real_2);
  auto result_2 = stan::math::ldexp(real_1, int_2);
  auto result_3 = stan::math::ldexp(int_1, real_2);
  auto result_4 = stan::math::ldexp(int_1, int_2);
  EXPECT_FLOAT_EQ(result_1, std::ldexp(real_1, real_2));
  EXPECT_FLOAT_EQ(result_2, std::ldexp(real_1, int_2));
  EXPECT_FLOAT_EQ(result_3, std::ldexp(int_1, real_2));
  EXPECT_FLOAT_EQ(result_4, std::ldexp(int_1, int_2));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, pow) {
  using namespace std::complex_literals;
  double real_1 = 0.4;
  double real_2 = 2.2;
  int int_1 = 2;
  int int_2 = 4;
  std::complex<double> complex_1 = 1. + 2.0i;
  std::complex<double> complex_2 = 2. + 3.0i;
  auto result_1 = stan::math::pow(real_1, real_2);
  auto result_2 = stan::math::pow(real_1, int_2);
  auto result_3 = stan::math::pow(real_1, complex_2);
  auto result_4 = stan::math::pow(int_1, real_2);
  auto result_5 = stan::math::pow(int_1, int_2);
  auto result_6 = stan::math::pow(int_1, complex_2);
  auto result_7 = stan::math::pow(complex_1, int_2);
  auto result_8 = stan::math::pow(complex_1, real_2);
  auto result_9 = stan::math::pow(complex_1, complex_2);
  EXPECT_FLOAT_EQ(result_1, std::pow(real_1, real_2));
  EXPECT_FLOAT_EQ(result_2, std::pow(real_1, int_2));
  EXPECT_EQ(result_3, std::pow(real_1, complex_2));
  EXPECT_FLOAT_EQ(result_4, std::pow(int_1, real_2));
  EXPECT_FLOAT_EQ(result_5, std::pow(int_1, int_2));
  EXPECT_EQ(result_6, std::pow(int_1, complex_2));
  EXPECT_EQ(result_7, std::pow(complex_1, int_2));
  EXPECT_EQ(result_8, std::pow(complex_1, real_2));
  EXPECT_EQ(result_9, std::pow(complex_1, complex_2));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}

TEST(PrimScalarSigTests, fma) {
  double real_1 = 0.4;
  double real_2 = 2.2;
  double real_3 = 4.2;
  int int_1 = 2;
  int int_2 = 4;
  int int_3 = 8;
  auto result_1 = stan::math::fma(real_1, real_2, real_3);
  auto result_2 = stan::math::fma(real_1, real_2, int_3);
  auto result_3 = stan::math::fma(real_1, int_2, real_3);
  auto result_4 = stan::math::fma(real_1, int_2, int_3);
  auto result_5 = stan::math::fma(int_1, real_2, real_3);
  auto result_6 = stan::math::fma(int_1, real_2, int_3);
  auto result_7 = stan::math::fma(int_1, int_2, real_3);
  auto result_8 = stan::math::fma(int_1, int_2, int_3);
  EXPECT_FLOAT_EQ(result_1, std::fma(real_1, real_2, real_3));
  EXPECT_FLOAT_EQ(result_2, std::fma(real_1, real_2, int_3));
  EXPECT_FLOAT_EQ(result_3, std::fma(real_1, int_2, real_3));
  EXPECT_FLOAT_EQ(result_4, std::fma(real_1, int_2, int_3));
  EXPECT_FLOAT_EQ(result_5, std::fma(int_1, real_2, real_3));
  EXPECT_FLOAT_EQ(result_6, std::fma(int_1, real_2, int_3));
  EXPECT_FLOAT_EQ(result_7, std::fma(int_1, int_2, real_3));
  EXPECT_FLOAT_EQ(result_8, std::fma(int_1, int_2, int_3));
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_stack_.size(), 0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_nochain_stack_.size(),
            0);
  EXPECT_EQ(stan::math::ChainableStack::instance_->var_alloc_stack_.size(), 0);
}
