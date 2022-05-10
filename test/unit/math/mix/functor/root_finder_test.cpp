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
  auto f = [](const auto& g, const auto& x) { return g * g * g * g * g - x; };
  auto f_g = [](const auto& g, const auto& x) { return 5 * g * g * g * g; };
  auto f_gg = [](const auto& g, const auto& x) { return 20 * g * g * g; };
  auto full_f = [&f, &f_g, &f_gg, guess, min, max](auto&& xx) {
    return root_finder(std::make_tuple(f, f_g, f_gg), guess, min, max, xx);
  };
  stan::test::expect_ad(full_f, x);
}
*/
TEST(MixFun, root_finder_beta_log) {
  auto beta_lpdf_deriv = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
    using stan::math::log;
    using stan::math::lgamma;
    using stan::math::log1m;
    auto val = -((beta * log1m(x) + (-2 + alpha) * log(x)
                  + log(1 - 2 * x - alpha + x * alpha + x * beta)
                  + lgamma(alpha + beta))
                 - (2 * log1m(x) + lgamma(alpha) + lgamma(beta)));
    return val;
  };
  auto f = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
    return stan::math::beta_lcdf(x, alpha, beta) - p;
  };
  auto f_g = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
    return stan::math::beta_lpdf(x, alpha, beta);
  };
  auto f_gg
      = [&beta_lpdf_deriv](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
          return beta_lpdf_deriv(x, alpha, beta, p);
        };
  double guess = .3;
  double min = 0.0000000001;
  double max = 1;
  auto full_f = [&f, &f_g, &f_gg, guess, min, max](auto&& alpha, auto&& beta,
                                                   auto&& p) {
    std::uintmax_t max_its = 1000;
    return stan::math::root_finder_hailey(std::make_tuple(f, f_g, f_gg), guess,
                                       min, max, alpha, beta, p);
  };
  double p = 0.45;
  double alpha = .5;
  double beta = .4;
  stan::test::expect_ad(full_f, alpha, beta, p);
}

TEST(MixFun, root_finder_beta) {
  auto beta_pdf_deriv = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
    using stan::math::lgamma;
    using stan::math::pow;
    auto val = -(pow((1 - x), beta) * pow(x, (-2 + alpha))
                 * (1. - 2. * x - alpha + x * alpha + x * beta)
                 * lgamma(alpha + beta))
               / (pow((-1 + x), 2) * lgamma(alpha) * lgamma(beta));
    return val;
  };
  auto beta_pdf = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
    using stan::math::lgamma;
    using stan::math::pow;
    auto val
        = (pow(x, alpha - 1) * pow((1 - x), beta - 1) * lgamma(alpha + beta))
          / ((lgamma(alpha) * lgamma(beta)));
    return val;
  };

  auto f = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
    return stan::math::beta_cdf(x, alpha, beta) - p;
  };
  double guess = .7;
  double min = 0.0000000001;
  double max = 1;
  auto full_f = [&f, &beta_pdf, &beta_pdf_deriv, guess, min, max](
                    auto&& alpha, auto&& beta, auto&& p) {
    std::uintmax_t max_its = 1000;
    return stan::math::root_finder_hailey(
        std::make_tuple(f, beta_pdf, beta_pdf_deriv), guess, min, max, alpha, beta, p);
  };
  double p = 0.8;
  double alpha = .4;
  double beta = .5;
  stan::test::expect_ad(full_f, alpha, beta, p);
}
