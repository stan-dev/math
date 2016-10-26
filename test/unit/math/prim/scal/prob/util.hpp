#ifndef TEST_UNIT_MATH_PRIM_SCAL_PROB_HPP
#define TEST_UNIT_MATH_PRIM_SCAL_PROB_HPP

#include <vector>
#include <boost/math/distributions.hpp>

/**
 * Uses a chi-squared test to assert that a vector of observed counts
 * is consistent with a vector of expected counts. Useful for testing RNGs.
*/
void assert_chi_squared(const std::vector<int>& counts, const std::vector<double>& expected, double tolerance) {
  int bins = counts.size();
	EXPECT_EQ(bins, expected.size());

  double chi = 0;
	for (int i=0; i<bins; ++i) {
    double discrepancy = expected[i] - counts[i];
		chi += discrepancy * discrepancy / expected[i];
  }

  boost::math::chi_squared dist(bins - 1);
	double chi_threshold = quantile(complement(dist, tolerance));

	EXPECT_TRUE(chi < chi_threshold);
}

/**
 * From a collection of samples and a list of quantiles, assumed ordered,
 * assert that the samples resemble draws from a distribution with those
 * quantiles, using a chi_squared goodness of fit test.
 */
void assert_matches_quantiles(const std::vector<double>& samples, const std::vector<double>& quantiles, double tolerance) {
  int N = samples.size();
	std::vector<double> mysamples = samples;
	std::sort(mysamples.begin(), mysamples.end());

	int K = quantiles.size();
	double expected_count = static_cast<double>(N) / K;

  std::vector<double> expected;
	for (int i=0; i<K; i++) {
	  expected.push_back(expected_count);
	}

	std::vector<int> counts(K);
  int current_index = 0;
	for (int i=0; i<N; ++i) {
    while (mysamples[i] >= quantiles[current_index]) {
      ++current_index;
    }
		++counts[current_index];
  }

	assert_chi_squared(counts, expected, tolerance);
}

#endif
