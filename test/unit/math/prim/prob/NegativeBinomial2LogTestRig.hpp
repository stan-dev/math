#ifndef TEST_UNIT_MATH_PRIM_PROB_NEGATIVE_BINOMIAL_2_LOG_TEST_RIG_HPP
#define TEST_UNIT_MATH_PRIM_PROB_NEGATIVE_BINOMIAL_2_LOG_TEST_RIG_HPP

#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/VectorRNGTestRig.hpp>
#include <cmath>

/*
 * NegativeBinomial2LogTestRig is a subclass of VectorRNGTestRig specialized for
 * testing random numbers generated from negative_binomial_2_log.
 *
 * Implements the function "pmf" and supplies values for test_points_
 *
 * The pmf callable wraps the probability mass function of the
 * negative_binomial_2_log.
 */
class NegativeBinomial2LogTestRig : public VectorIntRNGTestRig {
 public:
  NegativeBinomial2LogTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6},
                            {-0.5, 0.0, 0.1, 1.7}, {-2, 0, 1, 2}, {}, {},
                            {0.1, 1.1, 4.99}, {1, 2, 3}, {-3.0, -2.0, 0.0},
                            {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& eta, const T2& phi, const T3&,
                        T_rng& rng) const {
    return stan::math::neg_binomial_2_log_rng(eta, phi, rng);
  }

  template <typename T1>
  double pmf(int y, T1 eta, double phi, double) const {
    return std::exp(stan::math::neg_binomial_2_log_lpmf(y, eta, phi));
  }
};

#endif
