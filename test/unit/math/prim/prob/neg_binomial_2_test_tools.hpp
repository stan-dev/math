#ifndef STAN_TEST_UNIT_PRIM_PROB_NEG_BINOMIAL_2_TEST_TOOLS_HPP
#define STAN_TEST_UNIT_PRIM_PROB_NEG_BINOMIAL_2_TEST_TOOLS_HPP

#include <limits>
#include <cmath>
#include <stan/math/prim/prob/neg_binomial_2_lpmf.hpp>

namespace stan {
namespace test {
namespace neg_binomial_2_test_internal {
double phi_cutoff = stan::math::internal::neg_binomial_2_phi_cutoff;
double just_below_phi_cutoff = std::nextafter(phi_cutoff, 0);
double just_above_phi_cutoff
    = std::nextafter(phi_cutoff, std::numeric_limits<double>::infinity());
}  // namespace neg_binomial_2_test_internal
}  // namespace test
}  // namespace stan

#endif