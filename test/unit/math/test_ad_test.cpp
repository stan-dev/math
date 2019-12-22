#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>
#include <limits>
#include <type_traits>
#include <vector>

TEST(test_unit_math_test_ad, test_ad_unary) {
  Eigen::MatrixXd x(2, 2);
  x << 1.9, 0.3, 0.3, 1.7;

  auto g = [](const auto& u) { return stan::math::inverse(u); };
  stan::test::expect_ad(g, x);

  // Note:  the above is how tests need to be written;  the
  // functor is required because the following won't compile
  // because the template parameter F for expect_ad can't be deduced
  //   stan::test::expect_ad(inverse, x);
}

TEST(test_unit_math_test_ad, test_ad_binary) {
  auto g
      = [](const auto& u, const auto& v) { return stan::math::multiply(u, v); };

  Eigen::MatrixXd x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  double y = -3;
  stan::test::expect_ad(g, x, y);

  double y2 = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(g, x, y2);

  int y3 = -3;
  stan::test::expect_ad(g, x, y3);
}

template <typename T1, typename T2, typename T3>
stan::return_type_t<T1, T2, T3> arbitrary_ternary_function(T1 x1, T2 x2,
                                                           T3 x3) {
  stan::scalar_seq_view<T1> x1vec(x1);
  stan::scalar_seq_view<T2> x2vec(x2);
  stan::scalar_seq_view<T3> x3vec(x3);
  return x1vec[0] + x2vec[0] + x3vec[0];
}

TEST(test_unit_math_test_ad, test_ad_ternary) {
  auto g = [](const auto& x1, const auto& x2, const auto& x3) {
    return arbitrary_ternary_function(x1, x2, x3);
  };
  std::vector<double> v{1, 2, 3};
  double d = 1;
  int i = 1;
  stan::test::expect_ad(g, d, v, v);
  stan::test::expect_ad(g, i, v, v);
  stan::test::expect_ad(g, v, d, v);
  stan::test::expect_ad(g, v, i, v);
  stan::test::expect_ad(g, v, v, d);
  stan::test::expect_ad(g, v, v, i);
}

template <typename G, typename F, typename... Ts>
void expect_test_combinator_to_fail(const G& test_combinator, const F& f,
                                    const Ts&... xs) {
  // EXPECT_NOTFATAL_FAILURE succeeds if there is exactly one failure.
  // This function succeeds if there is one or more failures.
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(ADD_FAILURE() << "irrelevant failure message";
                              test_combinator(f, xs...), ""),
      "Expected: 1 non-fatal failure");
}

TEST(test_unit_math_test_ad, test_expect_test_combinator_to_fail) {
  auto combinator_that_succeeds = [](auto& f) {};
  auto combinator_that_fails = [](auto& f) { ADD_FAILURE() << ""; };
  auto combinator_that_fails_twice = [](auto& f) {
    ADD_FAILURE() << "";
    ADD_FAILURE() << "";
  };
  auto combinator_that_fails_fatally = [](auto& f) { FAIL() << ""; };
  auto combinator_that_fails_nonfatally_and_fatally = [](auto& f) {
    ADD_FAILURE() << "";
    FAIL() << "";
  };
  int fails = 0;
  auto combinator_that_tries_to_fail_fatally_twice = [&fails](auto& f) {
    ++fails;
    FAIL() << "";
    // this point not reached
    ++fails;
    FAIL() << "";
  };

  auto f = []() {};

  combinator_that_succeeds(f);
  EXPECT_NONFATAL_FAILURE(combinator_that_fails(f), "");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(combinator_that_fails_twice(f), ""),
      "Actual: 2 failures");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(combinator_that_fails_fatally(f), ""),
      "Fatal failure");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(combinator_that_fails_nonfatally_and_fatally(f),
                              ""),
      "Actual: 2 failures");

  EXPECT_NONFATAL_FAILURE(
      expect_test_combinator_to_fail(combinator_that_succeeds, f),
      "Actual: 0 failures");
  expect_test_combinator_to_fail(combinator_that_fails, f);
  expect_test_combinator_to_fail(combinator_that_fails_twice, f);
  expect_test_combinator_to_fail(combinator_that_fails_fatally, f);
  expect_test_combinator_to_fail(combinator_that_fails_nonfatally_and_fatally,
                                 f);
  fails = 0;
  expect_test_combinator_to_fail(combinator_that_tries_to_fail_fatally_twice,
                                 f);
  EXPECT_EQ(fails, 1);
}

template <typename F, typename... Ts>
void expect_expect_ad_failure(const F& f, const Ts&... xs) {
  auto g = [](auto& f, auto&... xs) {
    try {
      stan::test::expect_ad(f, xs...);
    } catch (...) {
      ADD_FAILURE();
    }
  };
  expect_test_combinator_to_fail(g, f, xs...);
}

// VECTORIZED UNARY FUNCTION THAT PASSES

// log10 is vectorized, so uses vectorized test
TEST(test_unit_math_test_ad, expect_ad_vectorized) {
  auto g = [](const auto& u) { return stan::math::log10(u); };

  stan::test::expect_ad_vectorized(g, 3.2);
}

// OVERLOAD THAT PASSES

// double overload matches template behavior and passes tests
template <typename T>
T f_match(const T& x) {
  return -2 * x;
}
double f_match(const double& x) { return -2 * x; }
TEST(test_unit_math_test_ad, match) {
  double x = 3.2;
  auto g = [](const auto& u) { return f_match(u); };
  stan::test::expect_ad(g, x);
}

// OVERLOAD THAT FAILS DUE TO MISMATCHED VALUE / DERIVATIVE

// double overload does not match template behavior for value
template <typename T>
T f_mismatch(const T& x) {
  return -2 * x;
}
stan::math::var f_mismatch(const stan::math::var& x) { return 2 * x; }
TEST(test_unit_math_test_ad, mismatch) {
  double x = 3.2;
  auto g = [](const auto& u) { return f_mismatch(u); };
  expect_expect_ad_failure(g, x);
}

// OVERLOAD THAT FAILS DUE TO MISMATCHED EXCEPTION CONDITIONS

// double overload does not match template behavior for exceptionsalue
template <typename T>
T f_misthrow(const T& x) {
  return -2 * x;
}
double f_misthrow(const double& x) {
  throw std::runtime_error("f_misthrow(double) called");
  return -2 * x;
}

TEST(test_unit_math_test_ad, misthrow) {
  double x = 1.73;
  auto h = [](const auto& u) { return f_misthrow(u); };
  expect_expect_ad_failure(h, x);
}

// ANOTHER OVERLOAD THAT FAILS DUE TO MISMATCHED EXCEPTION CONDITIONS

template <typename T1, typename T2>
stan::return_type_t<T1, T2> binary_function_misthrow(T1 x1, T2 x2) {
  stan::scalar_seq_view<T1> x1vec(x1);
  stan::scalar_seq_view<T2> x2vec(x2);
  return x1vec[0] + x2vec[0];
}

template <typename T1>
stan::return_type_t<T1> binary_function_misthrow(T1 x1, int x2) {
  throw std::runtime_error("binary_function_misthrow(T1, int) called");
  stan::scalar_seq_view<T1> x1vec(x1);
  stan::scalar_seq_view<int> x2vec(x2);
  return x1vec[0] + x2vec[0];
}

TEST(test_unit_math_test_ad, misthrow_binary) {
  auto g = [](const auto& x1, const auto& x2) {
    return binary_function_misthrow(x1, x2);
  };
  std::vector<double> v{1, 2, 3};
  double d = 1;
  int i = 1;
  stan::test::expect_ad(g, d, v);
  stan::test::expect_ad(g, i, v);
  stan::test::expect_ad(g, v, d);
  expect_expect_ad_failure(g, v, i);
}

// ANOTHER OVERLOAD THAT FAILS DUE TO MISMATCHED EXCEPTION CONDITIONS

template <typename T1, typename T2>
stan::return_type_t<T1, T2> binary_function_misthrow_2(T1 x1, T2 x2) {
  stan::scalar_seq_view<T1> x1vec(x1);
  stan::scalar_seq_view<T2> x2vec(x2);
  return x1vec[0] + x2vec[0];
}

template <typename T1>
stan::return_type_t<T1, stan::math::var> binary_function_misthrow_2(
    T1 x1, stan::math::var x2) {
  throw std::runtime_error("binary_function_misthrow(T1, var) called");
  stan::scalar_seq_view<T1> x1vec(x1);
  stan::scalar_seq_view<stan::math::var> x2vec(x2);
  return x1vec[0] + x2vec[0];
}

TEST(test_unit_math_test_ad, misthrow_binary_2) {
  auto g = [](const auto& x1, const auto& x2) {
    return binary_function_misthrow_2(x1, x2);
  };
  std::vector<double> v{1, 2, 3};
  double d = 1;
  stan::test::expect_ad(g, d, v);
  expect_expect_ad_failure(g, v, d);
}

// ANOTHER OVERLOAD THAT FAILS DUE TO MISMATCHED EXCEPTION CONDITIONS

template <typename T1, typename T2, typename T3>
stan::return_type_t<T1, T2, T3> ternary_function_misthrow(T1 x1, T2 x2, T3 x3) {
  stan::scalar_seq_view<T1> x1vec(x1);
  stan::scalar_seq_view<T2> x2vec(x2);
  stan::scalar_seq_view<T3> x3vec(x3);
  return x1vec[0] + x2vec[0] + x3vec[0];
}

template <typename T1, typename T2>
stan::return_type_t<T1, T2> ternary_function_misthrow(T1 x1, T2 x2, int x3) {
  throw std::runtime_error("ternary_function_misthrow(T1, T2, int) called");
  stan::scalar_seq_view<T1> x1vec(x1);
  stan::scalar_seq_view<T2> x2vec(x2);
  stan::scalar_seq_view<int> x3vec(x3);
  return x1vec[0] + x2vec[0] + x3vec[0];
}

TEST(test_unit_math_test_ad, misthrow_ternary) {
  auto g = [](const auto& x1, const auto& x2, const auto& x3) {
    return ternary_function_misthrow(x1, x2, x3);
  };
  std::vector<double> v{1, 2, 3};
  double d = 1;
  int i = 1;
  stan::test::expect_ad(g, d, v, v);
  stan::test::expect_ad(g, i, v, v);
  stan::test::expect_ad(g, v, d, v);
  stan::test::expect_ad(g, v, i, v);
  stan::test::expect_ad(g, v, v, d);
  expect_expect_ad_failure(g, v, v, i);
}

// ANOTHER OVERLOAD THAT FAILS DUE TO MISMATCHED EXCEPTION CONDITIONS

template <typename T1, typename T2, typename T3>
stan::return_type_t<T1, T2, T3> ternary_function_misthrow_2(T1 x1, T2 x2,
                                                            T3 x3) {
  stan::scalar_seq_view<T1> x1vec(x1);
  stan::scalar_seq_view<T2> x2vec(x2);
  stan::scalar_seq_view<T3> x3vec(x3);
  return x1vec[0] + x2vec[0] + x3vec[0];
}

template <typename T1, typename T2>
stan::return_type_t<T1, T2, stan::math::var> ternary_function_misthrow_2(
    T1 x1, T2 x2, stan::math::var x3) {
  throw std::runtime_error("ternary_function_misthrow(T1, T2, var) called");
  stan::scalar_seq_view<T1> x1vec(x1);
  stan::scalar_seq_view<T2> x2vec(x2);
  stan::scalar_seq_view<stan::math::var> x3vec(x3);
  return x1vec[0] + x2vec[0] + x3vec[0];
}

TEST(test_unit_math_test_ad, misthrow_ternary_2) {
  auto g = [](const auto& x1, const auto& x2, const auto& x3) {
    return ternary_function_misthrow_2(x1, x2, x3);
  };
  std::vector<double> v{1, 2, 3};
  double d = 1;
  stan::test::expect_ad(g, d, v, v);
  stan::test::expect_ad(g, v, d, v);
  expect_expect_ad_failure(g, v, v, d);
}

struct foo_fun {
  // must be static because operator() must be const
  static int calls_int_;
  static int calls_t_;

  template <typename T>
  T operator()(const T& x) const {
    ++foo_fun::calls_t_;
    return x / 2;
  }

  double operator()(const int& x) const {
    ++foo_fun::calls_int_;
    // return x / 2;  // fails if x is odd
    return x / 2.0;  // works for any x
  }
};

// need these or it won't link
int foo_fun::calls_int_ = -1;
int foo_fun::calls_t_ = -1;

TEST(test_unit_math_test_ad, integerGetsPassed) {
  // double arguments will not call int version
  foo_fun h;
  foo_fun::calls_int_ = 0;
  foo_fun::calls_t_ = 0;
  stan::test::expect_ad(h, 3.0);
  EXPECT_EQ(0, foo_fun::calls_int_);
  EXPECT_GT(foo_fun::calls_t_, 0);

  // int argument calls everything
  foo_fun g;
  foo_fun::calls_int_ = 0;
  foo_fun::calls_t_ = 0;
  stan::test::expect_ad(g, 3);
  EXPECT_GT(foo_fun::calls_int_, 0);
  EXPECT_GT(foo_fun::calls_t_, 0);

  foo_fun::calls_int_ = 0;
  foo_fun::calls_t_ = 0;
  stan::test::expect_common_unary(g);
  EXPECT_GT(foo_fun::calls_int_, 0);
  EXPECT_GT(foo_fun::calls_t_, 0);
}

struct bar_fun {
  // must be static because operator() must be const
  static int calls_t_;
  static int calls_int1_;
  static int calls_int2_;
  static int calls_int12_;

  static void reset() {
    bar_fun::calls_t_ = 0;
    bar_fun::calls_int1_ = 0;
    bar_fun::calls_int2_ = 0;
    bar_fun::calls_int12_ = 0;
  }

  double operator()(int x, int y) const {
    ++bar_fun::calls_int12_;
    return this->operator()(static_cast<double>(x), static_cast<double>(y));
  }
  double operator()(int x, double y) const {
    ++bar_fun::calls_int1_;
    return x * y;
  }
  template <typename T>
  T operator()(int x, const T& y) const {
    ++bar_fun::calls_int1_;
    return this->operator()(static_cast<double>(x), y);
  }

  double operator()(double x, int y) const {
    ++bar_fun::calls_int2_;
    return x * y;
  }
  double operator()(double x, double y) const { return x * y; }
  template <typename T>
  T operator()(double x, const T& y) const {
    return x * y;
  }

  template <typename T>
  T operator()(const T& x, int y) const {
    ++bar_fun::calls_int2_;
    return this->operator()(x, static_cast<double>(y));
  }
  template <typename T>
  T operator()(const T& x, double y) const {
    ++bar_fun::calls_t_;
    return x * y;
  }
  template <typename T1, typename T2>
  stan::promote_args_t<T1, T2> operator()(const T1& x, const T2& y) const {
    return x * y;
  }
};

// need these or it won't link
int bar_fun::calls_t_ = -1;
int bar_fun::calls_int1_ = -1;
int bar_fun::calls_int2_ = -1;
int bar_fun::calls_int12_ = -1;

TEST(test_unit_math_test_ad, integerGetsPassedBinary) {
  bar_fun f;

  bar_fun::reset();
  stan::test::expect_ad(f, 3.0, 2.7);
  EXPECT_GT(bar_fun::calls_t_, 0);
  EXPECT_EQ(bar_fun::calls_int1_, 0);
  EXPECT_EQ(bar_fun::calls_int2_, 0);
  EXPECT_EQ(bar_fun::calls_int12_, 0);

  bar_fun::reset();
  stan::test::expect_ad(f, 3, 2.7);
  EXPECT_GT(bar_fun::calls_t_, 0);
  EXPECT_GT(bar_fun::calls_int1_, 0);

  bar_fun::reset();
  stan::test::expect_ad(f, 2.7, 3);
  EXPECT_GT(bar_fun::calls_t_, 0);
  EXPECT_GT(bar_fun::calls_int2_, 0);

  bar_fun::reset();
  stan::test::expect_ad(f, 3, 4);
  EXPECT_GT(bar_fun::calls_t_, 0);
  EXPECT_GT(bar_fun::calls_int12_, 0);

  bar_fun::reset();
  stan::test::expect_common_binary(f);
  EXPECT_GT(bar_fun::calls_t_, 0);
  EXPECT_GT(bar_fun::calls_int1_, 0);
  EXPECT_GT(bar_fun::calls_int2_, 0);
  EXPECT_GT(bar_fun::calls_int12_, 0);
}

int baz_int = 0;
int baz_double = 0;
int baz_var = 0;
int baz_fvar = 0;

double baz(int x) {
  ++baz_int;
  return x / 2.0;
}
double baz(double x) {
  ++baz_double;
  return x / 2.0;
}
stan::math::var baz(const stan::math::var& x) {
  ++baz_var;
  return x / 2.0;
}
template <typename T>
stan::math::fvar<T> baz(const stan::math::fvar<T>& x) {
  ++baz_fvar;
  return x / 2.0;
}

struct baz_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return baz(x);
  }
};

template <typename T>
inline typename stan::math::apply_scalar_unary<baz_fun, T>::return_t baz(
    const T& x) {
  return stan::math::apply_scalar_unary<baz_fun, T>::apply(x);
}

TEST(test_unit_math_test_ad, integerGetsPassedVectorized) {
  auto h = [&](auto x) { return baz(x); };

  baz_int = 0;
  baz_double = 0;
  baz_var = 0;
  baz_fvar = 0;
  stan::test::expect_ad_vectorized(h, 0.2);
  EXPECT_EQ(0, baz_int);  // int doesn't get called
  EXPECT_GT(baz_double, 0);
  EXPECT_GT(baz_var, 0);
  EXPECT_GT(baz_fvar, 0);

  baz_int = 0;
  baz_double = 0;
  baz_var = 0;
  baz_fvar = 0;
  stan::test::expect_ad_vectorized(h, 1);
  EXPECT_GT(baz_int, 0);  // int version gets called
  EXPECT_GT(baz_double, 0);
  EXPECT_GT(baz_var, 0);
  EXPECT_GT(baz_fvar, 0);

  baz_int = 0;
  baz_double = 0;
  baz_var = 0;
  baz_fvar = 0;
  stan::test::expect_common_unary_vectorized(h);
  EXPECT_GT(baz_int, 0);  // int version gets called
  EXPECT_GT(baz_double, 0);
  EXPECT_GT(baz_var, 0);
  EXPECT_GT(baz_fvar, 0);

  baz_int = 0;
  baz_double = 0;
  baz_var = 0;
  baz_fvar = 0;
  stan::test::expect_common_nonzero_unary_vectorized(h);
  EXPECT_GT(baz_int, 0);  // int version gets called
  EXPECT_GT(baz_double, 0);
  EXPECT_GT(baz_var, 0);
  EXPECT_GT(baz_fvar, 0);
}

struct ternary_fun {
  // must be static because operator() must be const
  static int calls_int_;
  static int calls_int1_;
  static int calls_int2_;
  static int calls_int3_;
  static int calls_int12_;
  static int calls_int13_;
  static int calls_int23_;
  static int calls_int123_;

  static void reset() {
    ternary_fun::calls_int_ = 0;
    ternary_fun::calls_int1_ = 0;
    ternary_fun::calls_int2_ = 0;
    ternary_fun::calls_int3_ = 0;
    ternary_fun::calls_int12_ = 0;
    ternary_fun::calls_int13_ = 0;
    ternary_fun::calls_int23_ = 0;
    ternary_fun::calls_int123_ = 0;
  }

  double operator()(int x1, int x2, int x3) const {
    ++ternary_fun::calls_int123_;
    return this->operator()(static_cast<double>(x1), static_cast<double>(x2),
                            static_cast<double>(x3));
  }

  template <typename T1>
  T1 operator()(const T1& x1, int x2, int x3) const {
    ++ternary_fun::calls_int23_;
    return this->operator()(x1, static_cast<double>(x2),
                            static_cast<double>(x3));
  }

  template <typename T2>
  T2 operator()(int x1, const T2& x2, int x3) const {
    ++ternary_fun::calls_int13_;
    return this->operator()(static_cast<double>(x1), x2,
                            static_cast<double>(x3));
  }

  template <typename T3>
  T3 operator()(int x1, int x2, const T3& x3) const {
    ++ternary_fun::calls_int12_;
    return this->operator()(static_cast<double>(x1), static_cast<double>(x2),
                            x3);
  }

  template <typename T1, typename T2>
  stan::promote_args_t<T1, T2> operator()(const T1& x1, const T2& x2,
                                          int x3) const {
    ++ternary_fun::calls_int3_;
    return this->operator()(x1, x2, static_cast<double>(x3));
  }

  template <typename T1, typename T3>
  stan::promote_args_t<T1, T3> operator()(const T1& x1, int x2,
                                          const T3& x3) const {
    ++ternary_fun::calls_int2_;
    return this->operator()(x1, static_cast<double>(x2), x3);
  }

  template <typename T2, typename T3>
  stan::promote_args_t<T2, T3> operator()(int x1, const T2& x2,
                                          const T3& x3) const {
    ++ternary_fun::calls_int1_;
    return this->operator()(static_cast<double>(x1), x2, x3);
  }

  template <typename T1, typename T2, typename T3>
  stan::promote_args_t<T1, T2, T3> operator()(const T1& x1, const T2& x2,
                                              const T3& x3) const {
    ++ternary_fun::calls_int_;
    return x1 * x2 + x2 * x3 + x1 * x3;
  }
};

// need these or it won't link
int ternary_fun::calls_int_ = 0;
int ternary_fun::calls_int1_ = 0;
int ternary_fun::calls_int2_ = 0;
int ternary_fun::calls_int3_ = 0;
int ternary_fun::calls_int12_ = 0;
int ternary_fun::calls_int13_ = 0;
int ternary_fun::calls_int23_ = 0;
int ternary_fun::calls_int123_ = 0;

TEST(test_unit_math_test_ad, testAdTernaryIntPassed) {
  ternary_fun f;

  // { }
  ternary_fun::reset();
  stan::test::expect_ad(f, 1.0, 2.0, 3.0);
  EXPECT_LT(0, ternary_fun::calls_int_);
  EXPECT_EQ(0, ternary_fun::calls_int1_);
  EXPECT_EQ(0, ternary_fun::calls_int2_);
  EXPECT_EQ(0, ternary_fun::calls_int3_);
  EXPECT_EQ(0, ternary_fun::calls_int12_);
  EXPECT_EQ(0, ternary_fun::calls_int13_);
  EXPECT_EQ(0, ternary_fun::calls_int23_);
  EXPECT_EQ(0, ternary_fun::calls_int123_);

  // 1
  ternary_fun::reset();
  stan::test::expect_ad(f, 1, 2.0, 3.0);
  EXPECT_LT(0, ternary_fun::calls_int_);
  EXPECT_LT(0, ternary_fun::calls_int1_);
  EXPECT_EQ(0, ternary_fun::calls_int2_);
  EXPECT_EQ(0, ternary_fun::calls_int3_);
  EXPECT_EQ(0, ternary_fun::calls_int12_);
  EXPECT_EQ(0, ternary_fun::calls_int13_);
  EXPECT_EQ(0, ternary_fun::calls_int23_);
  EXPECT_EQ(0, ternary_fun::calls_int123_);

  // 2
  ternary_fun::reset();
  stan::test::expect_ad(f, 1.0, 2, 3.0);
  EXPECT_LT(0, ternary_fun::calls_int_);
  EXPECT_EQ(0, ternary_fun::calls_int1_);
  EXPECT_LT(0, ternary_fun::calls_int2_);
  EXPECT_EQ(0, ternary_fun::calls_int3_);
  EXPECT_EQ(0, ternary_fun::calls_int12_);
  EXPECT_EQ(0, ternary_fun::calls_int13_);
  EXPECT_EQ(0, ternary_fun::calls_int23_);
  EXPECT_EQ(0, ternary_fun::calls_int123_);

  // 3
  ternary_fun::reset();
  stan::test::expect_ad(f, 1.0, 2.0, 3);
  EXPECT_LT(0, ternary_fun::calls_int_);
  EXPECT_EQ(0, ternary_fun::calls_int1_);
  EXPECT_EQ(0, ternary_fun::calls_int2_);
  EXPECT_LT(0, ternary_fun::calls_int3_);
  EXPECT_EQ(0, ternary_fun::calls_int12_);
  EXPECT_EQ(0, ternary_fun::calls_int13_);
  EXPECT_EQ(0, ternary_fun::calls_int23_);
  EXPECT_EQ(0, ternary_fun::calls_int123_);

  // 1, 2
  ternary_fun::reset();
  stan::test::expect_ad(f, 1, 2, 3.0);
  EXPECT_LT(0, ternary_fun::calls_int_);
  EXPECT_LT(0, ternary_fun::calls_int1_);
  EXPECT_LT(0, ternary_fun::calls_int2_);
  EXPECT_EQ(0, ternary_fun::calls_int3_);
  EXPECT_LT(0, ternary_fun::calls_int12_);
  EXPECT_EQ(0, ternary_fun::calls_int13_);
  EXPECT_EQ(0, ternary_fun::calls_int23_);
  EXPECT_EQ(0, ternary_fun::calls_int123_);

  // 1, 3
  ternary_fun::reset();
  stan::test::expect_ad(f, 1, 2.0, 3);
  EXPECT_LT(0, ternary_fun::calls_int_);
  EXPECT_LT(0, ternary_fun::calls_int1_);
  EXPECT_EQ(0, ternary_fun::calls_int2_);
  EXPECT_LT(0, ternary_fun::calls_int3_);
  EXPECT_EQ(0, ternary_fun::calls_int12_);
  EXPECT_LT(0, ternary_fun::calls_int13_);
  EXPECT_EQ(0, ternary_fun::calls_int23_);
  EXPECT_EQ(0, ternary_fun::calls_int123_);

  // 2, 3
  ternary_fun::reset();
  stan::test::expect_ad(f, 1.0, 2, 3);
  EXPECT_LT(0, ternary_fun::calls_int_);
  EXPECT_EQ(0, ternary_fun::calls_int1_);
  EXPECT_LT(0, ternary_fun::calls_int2_);
  EXPECT_LT(0, ternary_fun::calls_int3_);
  EXPECT_EQ(0, ternary_fun::calls_int12_);
  EXPECT_EQ(0, ternary_fun::calls_int13_);
  EXPECT_LT(0, ternary_fun::calls_int23_);
  EXPECT_EQ(0, ternary_fun::calls_int123_);

  // 1, 2, 3
  ternary_fun::reset();
  stan::test::expect_ad(f, 1, 2, 3);
  EXPECT_LT(0, ternary_fun::calls_int_);
  EXPECT_LT(0, ternary_fun::calls_int1_);
  EXPECT_LT(0, ternary_fun::calls_int2_);
  EXPECT_LT(0, ternary_fun::calls_int3_);
  EXPECT_LT(0, ternary_fun::calls_int12_);
  EXPECT_LT(0, ternary_fun::calls_int13_);
  EXPECT_LT(0, ternary_fun::calls_int23_);
  EXPECT_LT(0, ternary_fun::calls_int123_);
}

template <typename T>
stan::return_type_t<T> prim_throws(const T& x) {
  throw std::runtime_error("prim_throws(T) called");
}

TEST(test_unit_math_test_ad, int_arguments_safe_when_prim_throws) {
  auto f = [](const auto& x) { return prim_throws(x); };
  stan::test::expect_ad(f, 1.0);
  // TODO(peterwicksstringfield) This is the problem with my first PR.
  expect_expect_ad_failure(f, 1);
}
