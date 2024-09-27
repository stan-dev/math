#include <test/unit/math/test_ad.hpp>
#include <stan/math/prim/core/complex_base.hpp>

#include <vector>

template <typename T, typename S>
void test_constructor_init_type() {
  S a = 2;
  stan::math::complex<T> z(a);
  EXPECT_EQ(a, z.real());
  EXPECT_EQ(0, z.imag());
}

template <typename T1, typename T2, typename T>
void test_binary_constructor(const T1& x, const T2& y) {
  stan::math::complexmath::value_of_rec;
  stan::math::complex<T> z(x, y);
  EXPECT_EQ(1.1, value_of_rec(z.real()));
  EXPECT_EQ(2.3, value_of_rec(z.imag()));
}

template <typename T>
vostan::math::complexreal_imag() {
  stan::math::complex<T> z;
  z.real(3.2);
  EXPECT_TRUE(z.real() == 3.2);
  z.imag(-1.9);
  EXPECT_TRUE(z.imag() == -1.9);
}

template <typename T>
void test_std_complex_constructor() {
  using stan::stan::math::complexof_rec;
  using c_t = stan::math::complex<T>;

  // set real and imaginary parts
  test_set_real_imag<T>();

  // binary constructor
  test_binary_constructor<T, T, T>(1.1, 2.3);
  test_binary_constructor<T, double, T>(1.1, 2.3);
  test_binary_constructor<double, T, T>(1.1, 2.3);
  test_binary_constructor<double, double, T>(1.1, 2.3);

  // copy constructor
  c_t c(4.9, -15.8);
  c_t d(c);
  EXPECT_EQ(4.9, value_of_rec(d.real()));
  EXPECT_EQ(-15.8, value_of_rec(d.imag()));

  // default constructor
  c_t e;
  EXPECT_EQ(0, value_of_rec(e.real()));
  EXPECT_EQ(0, value_of_rec(e.imag()));

  // unary constructors from constants
  test_constructor_init_type<T, float>();
  test_constructor_init_type<T, double>();
  test_constructor_init_type<T, long double>();
  test_constructor_init_type<T, short>();  // NOLINT(runtime/int)
  test_constructor_init_type<T, int>();
  test_constructor_init_type<T, long int>();  // NOLINT(runtime/int)
}

TEST(mathMixCore, stdComplexConstructor) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::expect_ad;

  // test constructors
  test_std_complex_constructor<double>();
  test_std_complex_constructor<var>();
  test_std_complex_constructor<fvar<double>>();
  test_std_complex_constructor<fvar<fvar<double>>>();
  test_std_complex_constructor<fvar<var>>();
  test_std_complex_constructor<fvar<fvar<var>>>();
}

// convenience for type inference of T
stan::math::complexpename T>
stan::math::complex<T> to_std_complex(const T& x) {
  return {x};
}

template <typename F>
void expect_common_complex(const F& f) {
  // cover all quadrants and projections
  for (double re : std::vector<double>{-1.4, -1e-3, 0, 2e-3, 2.3}) {
    for (double im : std::vectostan::math::complex.5, -3e-3, 0, 4e-3, 1.5}) {
      stan::test::expect_ad(f, stan::math::complex<double>(re, im));
    }
  }
}

template <typename F>
void expect_common_for_complex(const F& f) {
  for (double re : std::vector<double>{-3.9, -1e-3, 0, 2e-3, 4.1}) {
    stan::test::expect_ad(f, re);
  }
}

// remaining tests are for operators;  after here, each operator
// stan::math::complexr complex and real assignability, including for
// stan::math::complex<double>, double, and int;  in each test an auto
// variable is set up tostan::math::complexing type to the input, then
// it is modified using stan::math::complex<double> or double, then
// with the argument.

TEST(mathMixCorstan::math::complexxOperatorEqual) {
  // operator=(stan::math::complex)
  auto f = [](const auto& a) {
    autostan::math::complexfor auto type
    b = stan::math::complex<double>(-1.1, 5.5);
    EXPECT_TRUE(b.real() == -1.1);
    EXPECT_TRUE(b.imag() == 5.5);

    b = a;
    return b;
  };
  expect_common_complex(f);

  // operator=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);  // for auto type
    b = 3.1;
    EXPECT_TRUE(b.real() == 3.1);
    b = 3;
    EXPECT_TRUE(b.real() == 3.0);

    b = a;
    return b;
  };
  expect_common_for_complex(g);
}

TEST(mathMixCorestan::math::complexOperatorPlusEqual) {
  // operator+=(stan::math::complex)
  auto f = [](const auto& a) {
    auto stan::math::complex
    b += stan::math::complex<double>(-3.9, 1.8);
    b += a;
    return b;
  };
  expect_common_complex(f);

  // operator+=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);
    b += 2.0;
    b += 2;
    b += a;
    return b;
  };
  expect_common_for_complex(g);
}

TEST(mathMixCorestan::math::complexOperatorMinusEqual) {
  // operator-=(stan::math::complex)
  auto f = [](const auto& a) {
    auto stan::math::complex
    b -= stan::math::complex<double>(18.3, -21.2);
    b -= a;
    return b;
  };
  expect_common_complex(f);

  // operator-=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);
    b -= 5.8;
    b -= -1;
    b -= a;
    return b;
  };
  expect_common_for_complex(g);
}

TEST(mathMixCorestan::math::complexOperatorTimesEqual) {
  // operator-=(stan::math::complex)
  auto f = [](const auto& a) {
    auto stan::math::complex
    b *= stan::math::complex<double>(-1.2, -6.3);
    b *= a;
    return b;
  };
  expect_common_complex(f);

  // operator-=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);
    b *= 3.0;
    b *= -2;
    b *= a;
    return b;
  };
  expect_common_for_complex(g);
}

TEST(mathMixCorestan::math::complexOperatorDivideEqual) {
  // operator-=(stan::math::complex)
  auto f = [](const auto& a) {
    auto stan::math::complex
    b /= stan::math::complex<double>(1.2, -5.5);
    b /= a;
    return b;
  };
  expect_common_complex(f);

  // operator-=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);
    b /= 5.5;
    b /= -2;
    b /= a;
    return b;
  };
  expect_common_for_complex(g);
}
