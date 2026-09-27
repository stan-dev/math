#ifndef STAN_MATH_PRIM_FUN_INTERP1_CUBIC_HERMITE_HPP
#define STAN_MATH_PRIM_FUN_INTERP1_CUBIC_HERMITE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/meta.hpp>
#include <algorithm>
#include <iterator>
#include <tuple>

namespace stan {
namespace math {

/**
 * Precompute storage for cubic Hermite interpolation.
 *
 * This function sets up cubic Hermite interpolation by computing, for each
 * interval `[xk(i), xk(i + 1)]`, the coefficients of the cubic polynomial
 *   y(x) = a0 + a1 * dx + a2 * dx^2 + a3 * dx^3,
 * where `dx = x - xk(i)`. These coefficients are obtained by rewriting the
 * cubic Hermite basis representation
 *   y(t) = (1-t)^2 * (y0*(1+2*t) + s0*(x-x0)) + t^2 * (y1*(3-2*t) + dx*s1*(t-1))
 * (with `t = (x - x0) / (x1 - x0)`, `x0 = xk(i)`, `x1 = xk(i + 1)`,
 * `y0 = yk(i)`, `y1 = yk(i + 1)`, `s0 = dydxk(i)`, `s1 = dydxk(i + 1)`) as a
 * polynomial in `dx` rather than in `t`. See
 * boost/math/interpolators/detail/cubic_hermite_detail.hpp for background on
 * the cubic Hermite spline.
 *
 * @tparam EigVecX type of the knot locations
 * @tparam EigVecY type of the knot values
 * @param xk Strictly increasing interpolation knot locations.
 * @param yk Interpolation values at the knots.
 * @param dydxk Derivative values at the knots (user-provided).
 * @return A tuple containing a (K-1) x 4 coefficient matrix and the knot
 *   locations.
 * @throw std::invalid_argument if `xk`, `yk`, and `dydxk` do not have the same
 *   sizes
 * @throw std::domain_error if there are fewer than two knots, any input is
 *   not finite, or `xk` is not strictly increasing
 */
template <typename EigVecX, typename EigVecY,
          require_all_eigen_col_vector_t<EigVecX, EigVecY>* = nullptr>
inline std::tuple<
    Eigen::Matrix<return_type_t<EigVecX, EigVecY>, Eigen::Dynamic,
                  Eigen::Dynamic>,
    Eigen::Matrix<return_type_t<EigVecX>, Eigen::Dynamic, 1>>
interp1_cubic_hermite_setup(const EigVecX& xk, const EigVecY& yk,
                            const EigVecY& dydxk) {
  static constexpr const char* function = "interp1_cubic_hermite_setup";
  const auto& xk_ref = to_ref(xk);
  const auto& yk_ref = to_ref(yk);
  const auto& dydxk_ref = to_ref(dydxk);

  check_greater_or_equal(function, "number of knots", xk_ref.size(), 2);
  check_size_match(function, "xk", xk_ref.size(), "yk", yk_ref.size());
  check_size_match(function, "xk", xk_ref.size(), "dydxk", dydxk_ref.size());
  check_finite(function, "xk", xk_ref);
  check_finite(function, "yk", yk_ref);
  check_finite(function, "dydxk", dydxk_ref);
  check_ordered(function, "xk", xk_ref);

  using coef_t = return_type_t<EigVecX, EigVecY>;
  using xk_t = return_type_t<EigVecX>;
  const Eigen::Index K = xk_ref.size();
  Eigen::Matrix<coef_t, Eigen::Dynamic, Eigen::Dynamic> coef(K - 1, 4);

  for (Eigen::Index i = 0; i < K - 1; ++i) {
    const coef_t h = xk_ref.coeff(i + 1) - xk_ref.coeff(i);
    const coef_t y0 = yk_ref.coeff(i);
    const coef_t y1 = yk_ref.coeff(i + 1);
    const coef_t s0 = dydxk_ref.coeff(i);
    const coef_t s1 = dydxk_ref.coeff(i + 1);
    const coef_t h2 = h * h;
    const coef_t h3 = h2 * h;

    // Coefficients of y(x) = a0 + a1 * dx + a2 * dx^2 + a3 * dx^3, with
    // dx = x - xk(i), obtained by expressing the cubic Hermite basis
    // representation as a polynomial in dx rather than in
    // t = dx / h.
    coef(i, 0) = y0;
    coef(i, 1) = s0;
    coef(i, 2) = (3 * (y1 - y0) - h * (2 * s0 + s1)) / h2;
    coef(i, 3) = (2 * (y0 - y1) + h * (s0 + s1)) / h3;
  }

  Eigen::Matrix<xk_t, Eigen::Dynamic, 1> knots = xk_ref;
  return std::make_tuple(std::move(coef), std::move(knots));
}

/**
 * Evaluate cubic Hermite interpolation.
 *
 * This function evaluates the cubic Hermite interpolant built by
 * `interp1_cubic_hermite_setup` at the location `x`. It first locates the
 * interval `[xk(i), xk(i + 1)]` containing `x` (mirroring the search
 * performed in `cubic_hermite_detail::operator()` in
 * boost/math/interpolators/detail/cubic_hermite_detail.hpp), and then
 * evaluates the corresponding cubic polynomial
 *   y(x) = a0 + a1 * dx + a2 * dx^2 + a3 * dx^3,
 * with `dx = x - xk(i)`, using Horner's method. Because the coefficients
 * were derived so that this polynomial agrees exactly with the Hermite basis
 * representation on `[xk(i), xk(i + 1)]`, this reproduces the values that
 * would be returned by `cubic_hermite_detail::operator()`, including the
 * special handling at the right endpoint `x = xk(K - 1)`.
 *
 * @tparam T type of the evaluation location
 * @tparam Coef type of the coefficient matrix
 * @tparam Knots type of the knot locations
 * @param x Location at which to evaluate the interpolation.
 * @param W Tuple returned by `interp1_cubic_hermite_setup`.
 * @return The interpolated value at `x`.
 * @throw std::invalid_argument if `W` does not have the expected dimensions
 * @throw std::domain_error if `x` is not finite or is outside of the range
 *   spanned by the knots
 */
template <typename T, typename Coef, typename Knots,
          require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T, Coef, Knots> interp1_cubic_hermite_eval(
    const T& x, const std::tuple<Coef, Knots>& W) {
  static constexpr const char* function = "interp1_cubic_hermite_eval";
  const auto& coef = std::get<0>(W);
  const auto& xk = std::get<1>(W);

  check_size_match(function, "columns of coefficient matrix", coef.cols(), "4",
                   4);
  check_size_match(function, "rows of coefficient matrix plus one",
                   coef.rows() + 1, "size of xk", xk.size());
  check_finite(function, "x", x);

  const Eigen::Index K = xk.size();
  check_bounded(function, "x", x, xk.coeff(0), xk.coeff(K - 1));

  // Locate the interval i such that xk(i) <= x <= xk(i + 1). This mirrors
  // the std::upper_bound-based search in cubic_hermite_detail::operator():
  // `it` points to the first knot strictly greater than `x`, so `i` (one
  // less than the resulting index) is the interval containing `x`. When
  // `x` equals the last knot, `it == end`, giving `i == K - 1`; in that
  // case we fall back to the last interval, exactly as the reference
  // implementation special-cases `x == x_.back()`.
  const auto* begin = xk.data();
  const auto* end = begin + K;
  const auto it = std::upper_bound(begin, end, x);
  Eigen::Index i = std::distance(begin, it) - 1;
  if (i == K - 1) {
    --i;
  }

  const auto dx = x - xk.coeff(i);
  return coef(i, 0)
         + dx * (coef(i, 1) + dx * (coef(i, 2) + dx * coef(i, 3)));
}

}  // namespace math
}  // namespace stan

#endif
