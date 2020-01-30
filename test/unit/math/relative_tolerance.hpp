#ifndef TEST_UNIT_MATH_RELATIVE_TOLERANCE_HPP
#define TEST_UNIT_MATH_RELATIVE_TOLERANCE_HPP

#include <stan/math.hpp>
#include <algorithm>

namespace stan {
namespace test {

struct relative_tolerance {
  const double tol;
  const double tol_min;

  relative_tolerance(const double tol_)  // NOLINT
      : tol(tol_), tol_min(std::max(1e-14, tol_ * tol_)) {}

  relative_tolerance(const double tol_, const double tol_min_)
      : tol(tol_), tol_min(tol_min_) {}

  // Computes tolerance given an exact target value
  template <typename T1, require_stan_scalar_t<T1>...>
  double exact(const T1& x) const {
    using stan::math::fabs;
    return std::max(tol * fabs(x), tol_min);
  }

  // Computes tolerance given two inexact values
  template <typename T1, typename T2, require_all_stan_scalar_t<T1, T2>...>
  double inexact(const T1& x, const T1& y) const {
    using stan::math::fabs;
    return std::max(tol * 0.5 * (fabs(x) + fabs(y)), tol_min);
  }
};

const relative_tolerance default_relative_tolerance(1e-8, 1e-14);

}  // namespace test
}  // namespace stan
#endif
