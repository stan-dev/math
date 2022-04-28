#include <stan/math.hpp>
#include <stan/math/prim/functor/root_finder.hpp>
#include <gtest/gtest.h>


TEST(RevFunctor, root_finder) {
  using stan::math::root_finder;
  using stan::math::var;
  // return cube root of x using 1st and 2nd derivatives and Halley.
  //using namespace std;  // Help ADL of std functions.
  var x = 27;
  int exponent;
  std::frexp(x.val(), &exponent);                                // Get exponent of z (ignore mantissa).
  var guess = ldexp(1., exponent / 3);                    // Rough guess is to divide the exponent by three.
  var min = ldexp(0.5, exponent / 3);                     // Minimum possible value is half our guess.
  var max = ldexp(2., exponent / 3);                      // Maximum possible value is twice our guess.
  const int digits = std::numeric_limits<double>::digits;  // Maximum possible binary digits accuracy for type T.
  // digits used to control how accurate to try to make the result.
  int get_digits = static_cast<int>(digits * 0.4);    // Accuracy triples with each step, so stop when just
  // over one third of the digits are correct.
  std::uintmax_t maxit = 20;
  auto f = [](const auto& g, const auto x){ return g * g * g - x; };
  auto f_g = [](const auto& g, const auto x){ return 3 * g * g; };
  auto f_gg = [](const auto& g, const auto x){ return 6 * g; };
  var result = root_finder(f, f_g, f_gg, guess, min, max, x);
  result.grad();
  // val: 27 adj:0.037037
  std::cout << "\nval: " << x.val() << " adj:" << x.adj() << "\n";
}
