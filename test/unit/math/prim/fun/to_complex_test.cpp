#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

template <typename T, typename S>
void test_constructor(const T& re, const S& im) {
  std::complex<double> z1(re, im);
  std::complex<double> z2 = stan::math::to_complex(re, im);
  EXPECT_EQ(z1, z2);
}

TEST(mathPrimFunToComplex, isconstexpr) {
  using stan::math::to_complex;
  // using in static assert tests that constexpr really is const
  static_assert(std::real(to_complex(1, 2)) == 1, "Hello");
}

TEST(mathPrimFunToComplex, construction) {
  // test behavior for integer/double arg combos
  test_constructor(1, 2);
  test_constructor(1.0, 2);
  test_constructor(1, 2.0);
  test_constructor(1.0, 2.0);
}

TEST(MathFunctions, to_complex_vec) {
  Eigen::VectorXd in1(3);
  in1 << 1.2, 3.1, 0.8;
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;

  auto out1 = stan::math::to_complex(in1, in2);
  auto out2 = stan::math::to_complex(3, in2);

  Eigen::VectorXd three(3);
  three << 3.0, 3.0, 3.0;

  EXPECT_MATRIX_FLOAT_EQ(out1.real(), in1);
  EXPECT_MATRIX_FLOAT_EQ(out1.imag(), in2);
  EXPECT_MATRIX_FLOAT_EQ(out2.real(), three);
  EXPECT_MATRIX_FLOAT_EQ(out2.imag(), in2);
}
