#include <test/unit/math/test_ad.hpp>

struct CubedRootFinder {
  template <bool ReturnDeriv, typename T1, typename T2,
   std::enable_if_t<ReturnDeriv>* = nullptr>
  static auto run(T1&& g, T2&& x) {
    auto g_pow = g;
    auto second_deriv = 20 * g_pow;
    g_pow *= g;
    auto first_deriv = 3 * g_pow;
    g_pow *= g;
    auto func_val = g_pow - x;
    return std::make_tuple(func_val, first_deriv, second_deriv);
  }
  template <bool ReturnDeriv, typename T1, typename T2,
   std::enable_if_t<!ReturnDeriv>* = nullptr>
  static auto run(T1&& g, T2&& x) {
    return g * g * g - x;
  }
};

TEST(MixFun, root_finder_cubed) {
  using stan::math::root_finder;
  // return cube root of x using 1st and 2nd derivatives and Halley.
  // using namespace std;  // Help ADL of std functions.
  double x = 27;
  int exponent;
  // Get exponent of z (ignore mantissa).
  std::frexp(x, &exponent);
  // Rough guess is to divide the exponent by three.
  double guess = ldexp(1., exponent / 3);
  // Minimum possible value is half our guess.
  double min = ldexp(0.5, exponent / 3);
  // Maximum possible value is twice our guess.
  double max = ldexp(2., exponent / 3);
  // Maximum possible binary digits accuracy for type T.
  const int digits = std::numeric_limits<double>::digits;
  int get_digits = static_cast<int>(digits * 0.4);
  std::uintmax_t maxit = 20;
  auto full_f = [guess, min, max](auto&& xx) {
    return stan::math::root_finder_hailey<CubedRootFinder>(guess, min, max, xx);
  };
  stan::test::expect_ad(full_f, x);
}



struct FifthRootFinder {
  template <bool ReturnDeriv, typename T1, typename T2,
   std::enable_if_t<ReturnDeriv>* = nullptr>
  static auto run(T1&& g, T2&& x) {
    auto g_pow = g * g * g;
    auto second_deriv = 20 * g_pow;
    g_pow *= g;
    auto first_deriv = 5 * g_pow;
    g_pow *= g;
    auto func_val = g_pow - x;
    return std::make_tuple(func_val, first_deriv, second_deriv);
  }
  template <bool ReturnDeriv, typename T1, typename T2,
   std::enable_if_t<!ReturnDeriv>* = nullptr>
  static auto run(T1&& g, T2&& x) {
    return g * g * g * g * g - x;
  }
};

TEST(MixFun, root_finder_fifth) {
  using stan::math::root_finder;
  // return cube root of x using 1st and 2nd derivatives and Halley.
  // using namespace std;  // Help ADL of std functions.
  double x = 27;
  int exponent;
  // Get exponent of z (ignore mantissa).
  std::frexp(x, &exponent);
  // Rough guess is to divide the exponent by three.
  double guess = ldexp(1., exponent / 5);
  // Minimum possible value is half our guess.
  double min = ldexp(0.5, exponent / 5);
  // Maximum possible value is twice our guess.
  double max = ldexp(2., exponent / 5);
  // Maximum possible binary digits accuracy for type T.
  const int digits = std::numeric_limits<double>::digits;
  int get_digits = static_cast<int>(digits * 0.4);
  std::uintmax_t maxit = 20;
  auto f_hailey = [guess, min, max](auto&& x) {
    return stan::math::root_finder_hailey<FifthRootFinder>(guess, min, max, x);
  };
  stan::test::expect_ad(f_hailey, x);
}

struct BetaCdfRoot {
  template <bool ReturnDeriv, typename T1, typename T2, typename T3,
            typename T4, std::enable_if_t<ReturnDeriv>* = nullptr>
  static auto run(T1&& x, T2&& alpha, T3&& beta, T4&& p) {
    auto f_val = stan::math::inc_beta(alpha, beta, x) - p;
    auto beta_ab = stan::math::beta(alpha, beta);
    auto f_val_d = (stan::math::pow(1.0 - x,-1.0 + beta)*stan::math::pow(x,-1.0 + alpha)) / beta_ab;
    auto f_val_dd = (stan::math::pow(1 - x,-1 + beta)*stan::math::pow(x,-2 + alpha)*(-1 + alpha))/beta_ab -
   (stan::math::pow(1 - x,-2 + beta)*stan::math::pow(x,-1 + alpha)*(-1 + beta))/beta_ab;

    return std::make_tuple(f_val, f_val_d, f_val_dd);
  }
  template <bool ReturnDeriv, typename T1, typename T2, typename T3,
            typename T4, std::enable_if_t<!ReturnDeriv>* = nullptr>
  static auto run(T1&& x, T2&& alpha, T3&& beta, T4&& p) {
    return stan::math::inc_beta(alpha, beta, x) - p;
  }
};

TEST(MixFun, root_finder_beta) {
  constexpr double guess = .5;
  constexpr double min = 0;
  constexpr double max = 1;
  auto f_hailey = [](auto&& alpha, auto&& beta, auto&& p) {
    return stan::math::root_finder_hailey<BetaCdfRoot>(p, min, max, alpha, beta,
                                                       p);
  };
  double alpha = .5;
  double beta = .5;
  double p = .3;
  stan::test::expect_ad(f_hailey, alpha, beta, p);
}
