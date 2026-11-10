#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/incomplete_gamma_cdfs_reference.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace {

using stan::math::var;
using F3 = std::function<var(const var&, const var&, const var&)>;
using F2 = std::function<var(const var&, const var&)>;
using FP = std::function<var(int, const var&)>;

const std::map<std::string, F3>& three_arg() {
  namespace sm = stan::math;
  static const std::map<std::string, F3> f = {
      {"gamma_cdf",
       [](auto& y, auto& a, auto& b) { return sm::gamma_cdf(y, a, b); }},
      {"gamma_lcdf",
       [](auto& y, auto& a, auto& b) { return sm::gamma_lcdf(y, a, b); }},
      {"gamma_lccdf",
       [](auto& y, auto& a, auto& b) { return sm::gamma_lccdf(y, a, b); }},
      {"inv_gamma_cdf",
       [](auto& y, auto& a, auto& b) { return sm::inv_gamma_cdf(y, a, b); }},
      {"inv_gamma_lcdf",
       [](auto& y, auto& a, auto& b) { return sm::inv_gamma_lcdf(y, a, b); }},
      {"inv_gamma_lccdf",
       [](auto& y, auto& a, auto& b) { return sm::inv_gamma_lccdf(y, a, b); }},
      {"scaled_inv_chi_square_cdf",
       [](auto& y, auto& n, auto& s) {
         return sm::scaled_inv_chi_square_cdf(y, n, s);
       }},
      {"scaled_inv_chi_square_lcdf",
       [](auto& y, auto& n, auto& s) {
         return sm::scaled_inv_chi_square_lcdf(y, n, s);
       }},
      {"scaled_inv_chi_square_lccdf",
       [](auto& y, auto& n, auto& s) {
         return sm::scaled_inv_chi_square_lccdf(y, n, s);
       }},
  };
  return f;
}

const std::map<std::string, F2>& two_arg() {
  namespace sm = stan::math;
  static const std::map<std::string, F2> f = {
      {"chi_square_cdf",
       [](auto& y, auto& n) { return sm::chi_square_cdf(y, n); }},
      {"chi_square_lcdf",
       [](auto& y, auto& n) { return sm::chi_square_lcdf(y, n); }},
      {"chi_square_lccdf",
       [](auto& y, auto& n) { return sm::chi_square_lccdf(y, n); }},
      {"inv_chi_square_cdf",
       [](auto& y, auto& n) { return sm::inv_chi_square_cdf(y, n); }},
      {"inv_chi_square_lcdf",
       [](auto& y, auto& n) { return sm::inv_chi_square_lcdf(y, n); }},
      {"inv_chi_square_lccdf",
       [](auto& y, auto& n) { return sm::inv_chi_square_lccdf(y, n); }},
  };
  return f;
}

const std::map<std::string, FP>& poisson() {
  namespace sm = stan::math;
  static const std::map<std::string, FP> f = {
      {"poisson_cdf", [](int n, auto& l) { return sm::poisson_cdf(n, l); }},
      {"poisson_lcdf", [](int n, auto& l) { return sm::poisson_lcdf(n, l); }},
      {"poisson_lccdf", [](int n, auto& l) { return sm::poisson_lccdf(n, l); }},
  };
  return f;
}

/**
 * Relative error; a reference below 1e-300 (zero or subnormal, with few
 * significant bits) needs an absolute error below 1e-300.
 */
void expect_close(double expected, double actual, double tol,
                  const std::string& what) {
  ASSERT_TRUE(std::isfinite(actual)) << what << ": not finite, " << actual;
  if (std::fabs(expected) < 1e-300) {
    EXPECT_LE(std::fabs(actual - expected), 1e-300)
        << what << ": expected " << expected << ", got " << actual;
    return;
  }
  const double err = std::fabs(actual - expected) / std::fabs(expected);
  EXPECT_LT(err, tol) << what << ": relative error " << err << ", expected "
                      << expected << ", got " << actual;
}

}  // namespace

TEST_F(AgradRev, ProbDistributionsIncompleteGammaCdfs_reference) {
  namespace ref = stan::test::incomplete_gamma_cdfs;
  for (const auto& r : ref::rows()) {
    const std::string at = r.function + " at (" + std::to_string(r.y) + ", "
                           + std::to_string(r.p1) + ", " + std::to_string(r.p2)
                           + ")";
    std::vector<var> args;
    var lp;
    if (poisson().count(r.function)) {
      args = {var(r.p1)};
      lp = poisson().at(r.function)(static_cast<int>(r.y), args[0]);
    } else if (two_arg().count(r.function)) {
      args = {var(r.y), var(r.p1)};
      lp = two_arg().at(r.function)(args[0], args[1]);
    } else {
      ASSERT_TRUE(three_arg().count(r.function)) << r.function;
      args = {var(r.y), var(r.p1), var(r.p2)};
      lp = three_arg().at(r.function)(args[0], args[1], args[2]);
    }
    ASSERT_EQ(static_cast<int>(args.size()), r.n_grad) << at;
    lp.grad();
    expect_close(r.value, lp.val(), 1e-12, "value of " + at);
    for (int k = 0; k < r.n_grad; ++k) {
      expect_close(r.grad[k], args[k].adj(), 1e-9,
                   "gradient " + std::to_string(k) + " of " + at);
    }
    stan::math::recover_memory();
  }
}
