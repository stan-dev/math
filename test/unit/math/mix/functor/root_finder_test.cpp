#include <test/unit/math/test_ad.hpp>
/*
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
  auto f = [](const auto& g, const auto& x) { return g * g * g - x; };
  auto f_g = [](const auto& g, const auto& x) { return 3 * g * g; };
  auto f_gg = [](const auto& g, const auto& x) { return 6 * g; };
  auto full_f = [&f, &f_g, &f_gg, guess, min, max](auto&& xx) {
    return root_finder(std::make_tuple(f, f_g, f_gg), guess, min, max, xx);
  };
  stan::test::expect_ad(full_f, x);
}
*/


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

/*
struct BetaCdfRoot {
  template <bool ReturnDeriv, typename T1, typename T2, typename T3,
            typename T4, std::enable_if_t<ReturnDeriv>* = nullptr>
  static auto run(T1&& x, T2&& alpha, T3&& beta, T4&& p) {
    auto f = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
      return stan::math::beta_cdf(x, alpha, beta) - p;
    }(x, alpha, beta, p);
    auto beta_pdf_deriv = [](auto&& x, auto&& alpha, auto&& beta) {
      using stan::math::lgamma;
      using stan::math::pow;
      auto val = -(pow((1 - x), beta) * pow(x, (-2 + alpha))
                   * (1. - 2. * x - alpha + x * alpha + x * beta)
                   * lgamma(alpha + beta))
                 / (pow((-1 + x), 2) * lgamma(alpha) * lgamma(beta));
      return val;
    }(x, alpha, beta);
    auto beta_pdf = [](auto&& x, auto&& alpha, auto&& beta) {
      using stan::math::lgamma;
      using stan::math::pow;
      auto val
          = (pow(x, alpha - 1) * pow((1 - x), beta - 1) * lgamma(alpha + beta))
            / ((lgamma(alpha) * lgamma(beta)));
      return val;
    }(x, alpha, beta);
    return std::make_tuple(f, beta_pdf, beta_pdf_deriv);
  }
  template <bool ReturnDeriv, typename T1, typename T2, typename T3,
            typename T4, std::enable_if_t<!ReturnDeriv>* = nullptr>
  static auto run(T1&& x, T2&& alpha, T3&& beta, T4&& p) {
    return stan::math::beta_cdf(x, alpha, beta) - p;
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
*/
