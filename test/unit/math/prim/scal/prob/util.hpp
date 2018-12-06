#ifndef TEST_UNIT_MATH_PRIM_SCAL_PROB_HPP
#define TEST_UNIT_MATH_PRIM_SCAL_PROB_HPP

#include <boost/math/distributions.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

/**
 * Uses a chi-squared test to assert that a vector of observed counts
 * is consistent with a vector of expected counts. Useful for testing RNGs.
 */
void assert_chi_squared(const std::vector<int>& counts,
                        const std::vector<double>& expected, double tolerance) {
  int bins = counts.size();
  EXPECT_EQ(bins, expected.size());

  double chi = 0;
  for (int i = 0; i < bins; ++i) {
    double discrepancy = expected[i] - counts[i];
    chi += discrepancy * discrepancy / expected[i];
  }
  boost::math::chi_squared dist(bins - 1);
  double chi_threshold = quantile(complement(dist, tolerance));

  EXPECT_TRUE(chi < chi_threshold);
}

/**
 * From a collection of samples and a list of percentiles, assumed ordered,
 * assert that the samples resemble draws from a distribution with those
 * percentiles, using a chi_squared goodness of fit test.
 */
void assert_matches_cutpoints(const std::vector<double>& samples,
                              const std::vector<double>& pth_percentiles,
                              const std::vector<double>& ps,
                              double tolerance) {
  int N = samples.size();
  std::vector<double> mysamples = samples;
  std::sort(mysamples.begin(), mysamples.end());

  int K = pth_percentiles.size();
  assert(K == ps.size());
  std::vector<double> expected;
  for (int i = 0; i < K; i++) {
    assert(ps[i] >= 0 && ps[i] <= 1);
    expected.push_back(ps[i] * N);
  }

  std::vector<int> counts(K);
  size_t current_index = 0;
  for (int i = 0; i < N; ++i) {
    while (mysamples[i] >= pth_percentiles[current_index]) {
      ++current_index;
      EXPECT_TRUE(current_index < pth_percentiles.size());
    }
    ++counts[current_index];
  }
  assert_chi_squared(counts, expected, tolerance);
}

/**
 * From a collection of samples and a list of quantiles, assumed ordered,
 * assert that the samples resemble draws from a distribution with those
 * quantiles, using a chi_squared goodness of fit test.
 */
void assert_matches_quantiles(const std::vector<double>& samples,
                              const std::vector<double>& quantiles,
                              double tolerance) {
  int K = quantiles.size();
  std::vector<double> ps;
  for (int i = 0; i < K; ++i)
    ps.push_back(1.0 / K);

  assert_matches_cutpoints(samples, quantiles, ps, tolerance);
}
#endif
