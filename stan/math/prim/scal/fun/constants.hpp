#ifndef STAN_MATH_PRIM_SCAL_FUN_CONSTANTS_HPP
#define STAN_MATH_PRIM_SCAL_FUN_CONSTANTS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>
#include <boost/math/constants/constants.hpp>
#include <limits>

namespace stan {
namespace math {

using std::log;
using std::sqrt;

/**
 * The base of the natural logarithm,
 * \f$ e \f$.
 */
const double E = boost::math::constants::e<double>();

/**
 * The value of the square root of 2,
 * \f$ \sqrt{2} \f$.
 */
const double SQRT_2 = sqrt(2.0);

/**
 * The value of 1 over the square root of 2,
 * \f$ 1 / \sqrt{2} \f$.
 */
const double INV_SQRT_2 = inv(SQRT_2);

/**
 * The natural logarithm of 2,
 * \f$ \log 2 \f$.
 */
const double LOG_2 = log(2.0);

/**
 * The natural logarithm of 10,
 * \f$ \log 10 \f$.
 */
const double LOG_10 = log(10.0);

/**
 * Positive infinity.
 */
const double INFTY = std::numeric_limits<double>::infinity();

/**
 * Negative infinity.
 */
const double NEGATIVE_INFTY = -INFTY;

/**
 * (Quiet) not-a-number value.
 */
const double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();

/**
 * Smallest positive value.
 */
const double EPSILON = std::numeric_limits<double>::epsilon();

/**
 * Largest negative value (i.e., smallest absolute value).
 */
const double NEGATIVE_EPSILON = -EPSILON;

/**
 * Largest rate parameter allowed in Poisson RNG
 */
const double POISSON_MAX_RATE = pow(2.0, 30);

/**
 * Return the value of pi.
 *
 * @return Pi.
 */
inline double pi() { return boost::math::constants::pi<double>(); }

/**
 * Return the base of the natural logarithm.
 *
 * @return Base of natural logarithm.
 */
inline double e() { return E; }

/**
 * Return the square root of two.
 *
 * @return Square root of two.
 */
inline double sqrt2() { return SQRT_2; }

/**
 * Return natural logarithm of ten.
 *
 * @return Natural logarithm of ten.
 */
inline double log10() { return LOG_10; }

/**
 * Return positive infinity.
 *
 * @return Positive infinity.
 */
inline double positive_infinity() { return INFTY; }

/**
 * Return negative infinity.
 *
 * @return Negative infinity.
 */
inline double negative_infinity() { return NEGATIVE_INFTY; }

/**
 * Return (quiet) not-a-number.
 *
 * @return Quiet not-a-number.
 */
inline double not_a_number() { return NOT_A_NUMBER; }

/**
 * Returns the difference between 1.0 and the next value
 * representable.
 *
 * @return Minimum positive number.
 */
inline double machine_precision() { return EPSILON; }

const double SQRT_PI = sqrt(pi());

const double SQRT_2_TIMES_SQRT_PI = SQRT_2 * SQRT_PI;

const double TWO_OVER_SQRT_PI = 2.0 / SQRT_PI;

const double NEG_TWO_OVER_SQRT_PI = -TWO_OVER_SQRT_PI;

const double SQRT_TWO_PI = sqrt(2.0 * pi());

const double INV_SQRT_TWO_PI = inv(SQRT_TWO_PI);

const double LOG_PI = log(pi());

const double LOG_PI_OVER_FOUR = LOG_PI / 4.0;

const double LOG_SQRT_PI = log(SQRT_PI);

const double LOG_ZERO = log(0.0);

const double LOG_HALF = log(0.5);

const double NEG_LOG_TWO = -LOG_2;

const double NEG_LOG_SQRT_TWO_PI = -log(SQRT_TWO_PI);

const double NEG_LOG_PI = -LOG_PI;

const double NEG_LOG_SQRT_PI = -LOG_SQRT_PI;

const double NEG_LOG_TWO_OVER_TWO = -LOG_2 / 2.0;

const double LOG_TWO_PI = LOG_2 + LOG_PI;

const double NEG_LOG_TWO_PI = -LOG_TWO_PI;

const double LOG_EPSILON = log(EPSILON);
}  // namespace math
}  // namespace stan

#endif
