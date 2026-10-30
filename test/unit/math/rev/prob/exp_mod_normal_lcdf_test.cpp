#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/exp_mod_normal_lcdf_reference.hpp>
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <vector>

namespace {
const auto& ref_rows = stan::test::exp_mod_normal_lcdf_reference;

void expect_rel(double expected, double actual) {
  EXPECT_NEAR(expected, actual, 1e-12 * std::fabs(expected));
}

// Kind 0: lcdf, 1: lccdf, 2: cdf.
template <int Kind>
struct emg_fun {
  template <typename T_y, typename T_mu, typename T_sigma, typename T_lambda>
  static auto f(const T_y& y, const T_mu& mu, const T_sigma& sigma,
                const T_lambda& lambda) {
    if constexpr (Kind == 0) {
      return stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda);
    } else if constexpr (Kind == 1) {
      return stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda);
    } else {
      return stan::math::exp_mod_normal_cdf(y, mu, sigma, lambda);
    }
  }
};

// The reference value and the four partials of row r for Kind.
template <int Kind>
std::array<double, 5> expected(const std::array<double, 14>& r) {
  std::array<double, 5> e;
  const std::size_t k = Kind == 1 ? 9 : 4;
  for (std::size_t i = 0; i < 5; ++i) {
    e[i] = r[k + i];
  }
  if (Kind == 2) {
    e[0] = std::exp(r[4]);
    for (std::size_t i = 1; i < 5; ++i) {
      e[i] = e[0] * r[4 + i];
    }
  }
  return e;
}

template <int Kind>
void expect_rows_match_reference(
    const std::vector<std::array<double, 14>>& rows) {
  using stan::math::var;
  for (std::size_t n = 0; n < rows.size(); ++n) {
    SCOPED_TRACE(n);
    const auto& r = rows[n];
    const auto e = expected<Kind>(r);
    expect_rel(e[0], emg_fun<Kind>::f(r[0], r[1], r[2], r[3]));
    var y = r[0];
    var mu = r[1];
    var sigma = r[2];
    var lambda = r[3];
    var lp = emg_fun<Kind>::f(y, mu, sigma, lambda);
    lp.grad();
    expect_rel(e[0], lp.val());
    expect_rel(e[1], y.adj());
    expect_rel(e[2], mu.adj());
    expect_rel(e[3], sigma.adj());
    expect_rel(e[4], lambda.adj());
    stan::math::recover_memory();
  }
}

// All rows in one vectorized call of the lcdf or the lccdf: y and lambda
// vary, mu and sigma are shared scalars.
template <int Kind>
void expect_vectorized_sum_matches_reference() {
  using stan::math::var;
  std::vector<var> y;
  std::vector<var> lambda;
  double value = 0;
  double value_scale = 0;
  double d_mu = 0;
  double d_mu_scale = 0;
  double d_sigma = 0;
  double d_sigma_scale = 0;
  for (const auto& r : ref_rows) {
    const auto e = expected<Kind>(r);
    y.emplace_back(r[0]);
    lambda.emplace_back(r[3]);
    value += e[0];
    value_scale += std::fabs(e[0]);
    d_mu += e[2];
    d_mu_scale += std::fabs(e[2]);
    d_sigma += e[3];
    d_sigma_scale += std::fabs(e[3]);
  }
  var mu = ref_rows[0][1];
  var sigma = ref_rows[0][2];
  var lp = emg_fun<Kind>::f(y, mu, sigma, lambda);
  lp.grad();
  EXPECT_NEAR(value, lp.val(), 1e-12 * value_scale);
  EXPECT_NEAR(d_mu, mu.adj(), 1e-12 * d_mu_scale);
  EXPECT_NEAR(d_sigma, sigma.adj(), 1e-12 * d_sigma_scale);
  for (std::size_t n = 0; n < ref_rows.size(); ++n) {
    SCOPED_TRACE(n);
    const auto e = expected<Kind>(ref_rows[n]);
    expect_rel(e[1], y[n].adj());
    expect_rel(e[4], lambda[n].adj());
  }
}
}  // namespace

TEST_F(AgradRev, ProbDistributionsExpModNormal_lcdf_reference) {
  expect_rows_match_reference<0>(ref_rows);
}

TEST_F(AgradRev, ProbDistributionsExpModNormal_lccdf_reference) {
  expect_rows_match_reference<1>(ref_rows);
}

TEST_F(AgradRev, ProbDistributionsExpModNormal_cdf_reference) {
  expect_rows_match_reference<2>(ref_rows);
}

TEST_F(AgradRev, ProbDistributionsExpModNormal_lcdf_vectorized_reference) {
  expect_vectorized_sum_matches_reference<0>();
}

TEST_F(AgradRev, ProbDistributionsExpModNormal_lccdf_vectorized_reference) {
  expect_vectorized_sum_matches_reference<1>();
}

TEST_F(AgradRev, ProbDistributionsExpModNormal_deep_tail_reference) {
  const auto& rows = stan::test::exp_mod_normal_lcdf_deep_tail_reference;
  expect_rows_match_reference<0>(rows);
  expect_rows_match_reference<1>(rows);
  expect_rows_match_reference<2>(rows);
}

// Beyond z = -1.9e154, log F = -z^2 / 2 and its sigma partial overflow. The
// value is -inf (lcdf), 0 (lccdf) or 0 (cdf), never NaN; the other partials
// stay finite, and the cdf has zero partials.
TEST_F(AgradRev, ProbDistributionsExpModNormal_overflowing_tail) {
  using stan::math::var;
  for (double z : {-1e155, -1e200, -1e300}) {
    for (double v : {1e-300, 1.0, 1e10}) {
      SCOPED_TRACE(z);
      SCOPED_TRACE(v);
      {
        var y = z, mu = 0.0, sigma = 1.0, lambda = v;
        var lp = stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda);
        lp.grad();
        EXPECT_EQ(lp.val(), stan::math::NEGATIVE_INFTY);
        // d log F / dz = -z + v + O(1 / z) here
        EXPECT_FLOAT_EQ(y.adj(), -z);
        EXPECT_FLOAT_EQ(mu.adj(), z);
        EXPECT_TRUE(std::isfinite(lambda.adj()));
        stan::math::recover_memory();
      }
      {
        var y = z, mu = 0.0, sigma = 1.0, lambda = v;
        var lp = stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda);
        lp.grad();
        EXPECT_EQ(lp.val(), 0.0);
        EXPECT_TRUE(std::isfinite(y.adj()));
        EXPECT_TRUE(std::isfinite(sigma.adj()));
        EXPECT_TRUE(std::isfinite(lambda.adj()));
        stan::math::recover_memory();
      }
      {
        var y = z, mu = 0.0, sigma = 1.0, lambda = v;
        var lp = stan::math::exp_mod_normal_cdf(y, mu, sigma, lambda);
        lp.grad();
        EXPECT_EQ(lp.val(), 0.0);
        EXPECT_EQ(y.adj(), 0.0);
        EXPECT_EQ(mu.adj(), 0.0);
        EXPECT_EQ(sigma.adj(), 0.0);
        EXPECT_EQ(lambda.adj(), 0.0);
        stan::math::recover_memory();
      }
    }
  }
}
