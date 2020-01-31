#ifndef TEST_UNIT_MATH_RELATIVE_TOLERANCE_HPP
#define TEST_UNIT_MATH_RELATIVE_TOLERANCE_HPP

#include <stan/math.hpp>
#include <algorithm>

namespace stan {
namespace test {

/**
 * Struct holding information about relative tolerance and the minimal absolute
 * tolerance that should be tested against. The final tolerance is computed as
 * max(tol * fabs(x), tol_min)  
 * Where x is either an exact value to be tested against or average of 
 * two inexact values to be compared.
 */
struct relative_tolerance {
  /** 
   * The relative tolerance
   */
  double tol;
  /** 
   * The minimal absolute tolerance
   */
  double tol_min;

  /**
   * Construct with default tolerances
   */
  relative_tolerance() : tol(1e-8), tol_min(1e-14) {}

  /**
   * Construct with default tol_min (max(tol_ * tol_, 1e-14))
   * @param tol_ the relative tolerance
   */
  relative_tolerance(const double tol_)  // NOLINT
      : tol(tol_), tol_min(std::max(tol_ * tol_, 1e-14)) {}

  /**
   * Construct fully specified
   * @param[in] tol_ the relative tolerance
   * @param[in] tol_min_ the minimum absolute tolerance
   */
  relative_tolerance(const double tol_, const double tol_min_)
      : tol(tol_), tol_min(tol_min_) {}

  /**
   * Computes tolerance around an exact target value.
   * 
   * @tparam T1 the type of argument, must be scalar
   * @param[in] x the target value
   * @return tolerance
   */
  template <typename T1, require_stan_scalar_t<T1>...>
  double exact(const T1& x) const {
    using stan::math::fabs;
    return std::max(tol * fabs(x), tol_min);
  }

  /**
   * Computes tolerance for comparing two inexact values.
   * 
   * @tparam T1 the type of first argument, must be scalar
   * @tparam T2 the type of second argument, must be scalar
   * @param[in] x first value that will be compared
   * @param[in] y second value that will be compared
   * @return tolerance
   */
  template <typename T1, typename T2, require_all_stan_scalar_t<T1, T2>...>
  double inexact(const T1& x, const T2& y) const {
    using stan::math::fabs;
    return std::max(tol * 0.5 * (fabs(x) + fabs(y)), tol_min);
  }
};

}  // namespace test
}  // namespace stan
#endif
