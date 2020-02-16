#include <test/unit/math/test_ad.hpp>
#include <complex>

template <typename T, typename S>
void test_constructor_init_type() {
  S a = 2;
  std::complex<T> z(a);
  EXPECT_EQ(a, z.real());
  EXPECT_EQ(0, z.imag());
}

template <typename T1, typename T2, typename T>
void test_binary_constructor(const T1& x, const T2& y) {
  using stan::math::value_of_rec;
  std::complex<T> z(x, y);
  EXPECT_EQ(1.1, value_of_rec(z.real()));
  EXPECT_EQ(2.3, value_of_rec(z.imag()));
}

template <typename T>
void test_std_complex_constructor() {
  using stan::math::value_of_rec;
  using c_t = std::complex<T>;

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
  test_constructor_init_type<T, short>();
  test_constructor_init_type<T, int>();
  test_constructor_init_type<T, long int>();
}

template <typename T>
std::vector<T> my_to_vector(const T& x) {
  return {x};
}
template <typename T>
std::vector<T> my_to_vector(const std::vector<T>& x) {
  return x;
}
template <typename T>
std::vector<T> my_to_vector(const std::complex<T>& z) {
  return {z.real(), z.imag()};
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

template <typename T>
std::complex<T> to_std_complex(const T& x) {
  return {x};
}

TEST(mathMixCore, stdComplexOperatorEqual) {
  using stan::test::expect_ad;

  // operator=(std::complex)
  auto f = [](const auto& a) {
    auto b = a;  // for auto type
    b = a;
    return b;
  };
  expect_ad(f, std::complex<double>(2.3, 1.9));

  // operator=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);  // for auto type
    b = a;
    return b;
  };
  expect_ad(g, 2.3);
}

TEST(mathMixCore, stdComplexOperatorPlusEqual) {
  using stan::test::expect_ad;

  // operator+=(std::complex)
  auto f = [](const auto& a) {
    auto b = a;
    b += a;
    return b;
  };
  expect_ad(f, std::complex<double>(2.3, 1.9));

  // operator+=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);
    b += a;
    return b;
  };
  expect_ad(g, 2.3);
}

TEST(mathMixCore, stdComplexOperatorMinusEqual) {
  using stan::test::expect_ad;

  // operator-=(std::complex)
  auto f = [](const auto& a) {
    auto b = a;
    b -= a;
    return b;
  };
  expect_ad(f, std::complex<double>(2.3, 1.9));

  // operator-=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);
    b -= a;
    return b;
  };
  expect_ad(g, 2.3);
}

TEST(mathMixCore, stdComplexOperatorTimesEqual) {
  using stan::test::expect_ad;

  // operator-=(std::complex)
  auto f = [](const auto& a) {
    auto b = a;
    b *= a;
    return b;
  };
  expect_ad(f, std::complex<double>(2.3, 1.9));

  // operator-=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);
    b *= a;
    return b;
  };
  expect_ad(g, 2.3);
}

TEST(mathMixCore, stdComplexOperatorDivideEqual) {
  using stan::test::expect_ad;

  // operator-=(std::complex)
  auto f = [](const auto& a) {
    auto b = a;
    b /= a;
    return b;
  };
  expect_ad(f, std::complex<double>(2.3, 1.9));

  // operator-=(Arith)
  auto g = [](const auto& a) {
    auto b = to_std_complex(a);
    b /= a;
    return b;
  };
  expect_ad(g, 2.3);
}
