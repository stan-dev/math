#ifndef STAN_MATH_PRIM_FUN_CONSTANTS_HPP
#define STAN_MATH_PRIM_FUN_CONSTANTS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

// TODO(anyone) Use constexpr when moving to C++17

/**
 * Return the base of the natural logarithm.
 *
 * @return Base of natural logarithm.
 */
static constexpr double e() { return boost::math::constants::e<double>(); }

/**
 * Return the Euler's gamma constant.
 *
 * @return Euler's Gamma.
 */
static constexpr double egamma() {
  return boost::math::constants::euler<double>();
}

/**
 * Return the value of pi.
 *
 * @return Pi.
 */
static constexpr double pi() { return boost::math::constants::pi<double>(); }

/**
 * Smallest positive value.
 */
static constexpr double EPSILON = std::numeric_limits<double>::epsilon();

/**
 * Positive infinity.
 */
static constexpr double INFTY = std::numeric_limits<double>::infinity();

/**
 * Negative infinity.
 */
static constexpr double NEGATIVE_INFTY = -INFTY;

/**
 * (Quiet) not-a-number value.
 */
static constexpr double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();

/**
 * Twice the value of \f$ \pi \f$,
 * \f$ 2\pi \f$.
 */
static constexpr double TWO_PI = boost::math::constants::two_pi<double>();

/**
 * The natural logarithm of 0,
 * \f$ \log 0 \f$.
 */
static constexpr double LOG_ZERO = -INFTY;

/**
 * The natural logarithm of machine precision \f$ \epsilon \f$,
 * \f$ \log \epsilon \f$.
 */
const double LOG_EPSILON = std::log(EPSILON);

/**
 * The natural logarithm of 2,
 * \f$ \log 2 \f$.
 */
static constexpr double LOG_TWO = boost::math::constants::ln_two<double>();

/**
 * The natural logarithm of \f$ \pi \f$,
 * \f$ \log \pi \f$.
 */
static constexpr double LOG_PI = 1.14472988584940017414342735135;

/**
 * The natural logarithm of 0.5,
 * \f$ \log 0.5 = \log 1 - \log 2 \f$.
 */
static constexpr double LOG_HALF = -LOG_TWO;

/**
 * The natural logarithm of 2 plus the natural logarithm of \f$ \pi \f$,
 * \f$ \log(2\pi) \f$.
 */
static constexpr double LOG_TWO_PI = LOG_TWO + LOG_PI;

/**
 * The value of one quarter the natural logarithm of \f$ \pi \f$,
 * \f$ \log(\pi) / 4 \f$.
 */
static constexpr double LOG_PI_OVER_FOUR = 0.25 * LOG_PI;

/**
 * The natural logarithm of the square root of \f$ \pi \f$,
 * \f$ \log(sqrt{\pi}) \f$.
 */
static constexpr double LOG_SQRT_PI = LOG_PI / 2;

/**
 * The natural logarithm of 10,
 * \f$ \log 10 \f$.
 */
static constexpr double LOG_TEN = boost::math::constants::ln_ten<double>();

/**
 * The value of the square root of 2,
 * \f$ \sqrt{2} \f$.
 */
static constexpr double SQRT_TWO = boost::math::constants::root_two<double>();

/**
 * The value of the square root of \f$ \pi \f$,
 * \f$ \sqrt{\pi} \f$.
 */
static constexpr double SQRT_PI = boost::math::constants::root_pi<double>();

/**
 * The value of the square root of \f$ 2\pi \f$,
 * \f$ \sqrt{2\pi} \f$.
 */
static constexpr double SQRT_TWO_PI
    = boost::math::constants::root_two_pi<double>();

/**
 * The square root of 2 divided by the square root of \f$ \pi \f$,
 * \f$ \sqrt{2} / \sqrt{\pi} \f$.
 */
static constexpr double SQRT_TWO_OVER_SQRT_PI = SQRT_TWO / SQRT_PI;

/**
 * The value of 1 over the square root of 2,
 * \f$ 1 / \sqrt{2} \f$.
 */
static constexpr double INV_SQRT_TWO
    = boost::math::constants::one_div_root_two<double>();

/**
 * The value of 1 over the square root of \f$ \pi \f$,
 * \f$ 1 / \sqrt{\pi} \f$.
 */
static constexpr double INV_SQRT_PI
    = boost::math::constants::one_div_root_pi<double>();

/**
 * The value of 1 over the square root of \f$ 2\pi \f$,
 * \f$ 1 / \sqrt{2\pi} \f$.
 */
static constexpr double INV_SQRT_TWO_PI
    = boost::math::constants::one_div_root_two_pi<double>();

/**
 * The value of 2 over the square root of \f$ \pi \f$,
 * \f$ 2 / \sqrt{\pi} \f$.
 */
static constexpr double TWO_OVER_SQRT_PI
    = boost::math::constants::two_div_root_pi<double>();

/**
 * The value of half the natural logarithm 2,
 * \f$ \log(2) / 2 \f$.
 */
static constexpr double HALF_LOG_TWO = 0.5 * LOG_TWO;

/**
 * The value of half the natural logarithm \f$ 2\pi \f$,
 * \f$ \log(2\pi) / 2 \f$.
 */
static constexpr double HALF_LOG_TWO_PI = 0.5 * LOG_TWO_PI;

/**
 * The value of minus the natural logarithm of the square root of \f$ 2\pi \f$,
 * \f$ -\log(\sqrt{2\pi}) \f$.
 */
const double NEG_LOG_SQRT_TWO_PI = -std::log(SQRT_TWO_PI);

/**
 * Largest rate parameter allowed in Poisson RNG
 */
const double POISSON_MAX_RATE = std::pow(2.0, 30);

/**
 * Return positive infinity.
 *
 * @return Positive infinity.
 */
static constexpr inline double positive_infinity() { return INFTY; }

/**
 * Return negative infinity.
 *
 * @return Negative infinity.
 */
static constexpr inline double negative_infinity() { return NEGATIVE_INFTY; }

/**
 * Return (quiet) not-a-number.
 *
 * @return Quiet not-a-number.
 */
static constexpr inline double not_a_number() { return NOT_A_NUMBER; }

/**
 * Returns the difference between 1.0 and the next value
 * representable.
 *
 * @return Minimum positive number.
 */
static constexpr inline double machine_precision() { return EPSILON; }

/**
 * Returns the natural logarithm of ten.
 *
 * @return Natural logarithm of ten.
 */
static constexpr inline double log10() { return LOG_TEN; }

/**
 * Returns the square root of two.
 *
 * @return Square root of two.
 */
static constexpr inline double sqrt2() { return SQRT_TWO; }

}  // namespace math
}  // namespace stan

#endif
