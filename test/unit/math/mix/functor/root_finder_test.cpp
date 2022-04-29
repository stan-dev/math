#include <test/unit/math/test_ad.hpp>

TEST(MixFun, root_finder_cubed) {
  using stan::math::root_finder;
  // return cube root of x using 1st and 2nd derivatives and Halley.
  //using namespace std;  // Help ADL of std functions.
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
  /* digits used to control how accurate to try to make the result.
   * Accuracy triples with each step, so stop when just
   * over one third of the digits are correct.
   */
  int get_digits = static_cast<int>(digits * 0.4);
  std::uintmax_t maxit = 20;
  auto f = [](const auto& g, const auto& x){ return g * g * g - x; };
  auto f_g = [](const auto& g, const auto& x){ return 3 * g * g; };
  auto f_gg = [](const auto& g, const auto& x){ return 6 * g; };
  auto full_f = [&f, &f_g, &f_gg, guess, min, max](auto&& xx) {
    using x_t = std::decay_t<decltype(xx)>;
    return root_finder(std::make_tuple(f, f_g, f_gg),
     x_t(guess),
     x_t(min),
     x_t(max), xx);
  };
  stan::test::expect_ad(full_f, x);
}

TEST(MixFun, root_finder_fifth) {
  using stan::math::root_finder;
  // return cube root of x using 1st and 2nd derivatives and Halley.
  //using namespace std;  // Help ADL of std functions.
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
  /* digits used to control how accurate to try to make the result.
   * Accuracy triples with each step, so stop when just
   * over one third of the digits are correct.
   */
  int get_digits = static_cast<int>(digits * 0.4);
  std::uintmax_t maxit = 20;
  auto f = [](const auto& g, const auto& x){ return g * g * g * g * g - x; };
  auto f_g = [](const auto& g, const auto& x){ return 5 * g * g * g * g; };
  auto f_gg = [](const auto& g, const auto& x){ return 20 * g * g * g; };
  auto full_f = [&f, &f_g, &f_gg, guess, min, max](auto&& xx) {
    using x_t = std::decay_t<decltype(xx)>;
    return root_finder(std::make_tuple(f, f_g, f_gg),
     x_t(guess),
     x_t(min),
     x_t(max), xx);
  };
  stan::test::expect_ad(full_f, x);
}
