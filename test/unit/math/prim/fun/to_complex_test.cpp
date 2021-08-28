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
