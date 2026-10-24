#ifndef STAN_MATH_PRIM_FUN_INTERPOLATE_1D_CUBIC_HERMITE_HPP
#define STAN_MATH_PRIM_FUN_INTERPOLATE_1D_CUBIC_HERMITE_HPP

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
 * This function validates the input knot locations, values, and derivatives,
 * then stores them in a tuple for use with `interpolate_1d_cubic_spline_eval`.
 *
 * @tparam EigVecX type of the knot locations
 * @tparam EigVecY type of the knot values
 * @tparam EigVecDyDx type of the derivative values
 * @param xk Strictly increasing interpolation knot locations.
 * @param yk Interpolation values at the knots.
 * @param dydxk Derivative values at the knots (user-provided).
 * @return A tuple containing:
 *   - Knot locations `xk`
 *   - Knot values `yk`
 *   - Derivative values `dydxk`
 * @throw std::invalid_argument if `xk`, `yk`, and `dydxk` do not have the same
 *   sizes
 * @throw std::domain_error if there are fewer than two knots, any input is
 *   not finite, or `xk` is not strictly increasing
 */
template <typename EigVecX, typename EigVecY, typename EigVecDyDx,
          require_all_eigen_col_vector_t<EigVecX, EigVecY, EigVecDyDx>* = nullptr>
inline std::tuple<
    Eigen::Matrix<return_type_t<EigVecX>, Eigen::Dynamic, 1>,
    Eigen::Matrix<return_type_t<EigVecY>, Eigen::Dynamic, 1>,
    Eigen::Matrix<return_type_t<EigVecDyDx>, Eigen::Dynamic, 1>>
interpolate_1d_cubic_hermite_setup(const EigVecX& xk, const EigVecY& yk,
                            const EigVecDyDx& dydxk) {
  static constexpr const char* function = "interpolate_1d_cubic_hermite_setup";
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

  using xk_t = return_type_t<EigVecX>;
  using yk_t = return_type_t<EigVecY>;
  using dydxk_t = return_type_t<EigVecDyDx>;

  Eigen::Matrix<xk_t, Eigen::Dynamic, 1> xk_out = xk_ref;
  Eigen::Matrix<yk_t, Eigen::Dynamic, 1> yk_out = yk_ref;
  Eigen::Matrix<dydxk_t, Eigen::Dynamic, 1> dydxk_out = dydxk_ref;

  return std::make_tuple(std::move(xk_out), std::move(yk_out), std::move(dydxk_out));
}
/**
 * Evaluate cubic Hermite interpolation.
 *
 * This function evaluates the cubic Hermite interpolant at location `x` using
 * the standard scaled cubic Hermite basis functions:
 *   y(x) = (1-t)^2 * (y0*(1+2*t) + h*s0*t) + t^2 * (y1*(3-2*t) + h*s1*(t-1))
 * where `t = (x - xk(i)) / h`, `h = xk(i+1) - xk(i)`, `y0 = yk(i)`,
 * `y1 = yk(i+1)`, `s0 = dydxk(i)`, `s1 = dydxk(i+1)`.
 *
 * @tparam T type of the evaluation location
 * @tparam XkType type of the knot locations
 * @tparam YkType type of the knot values
 * @tparam DydxkType type of the derivative values
 * @param x Location at which to evaluate the interpolation.
 * @param W Tuple returned by `interpolate_1d_cubic_hermite_setup` containing (xk, yk, dydxk).
 * @return The interpolated value at `x`.
 * @throw std::invalid_argument if `W` does not have the expected structure
 * @throw std::domain_error if `x` is not finite or is outside of the range
 *   spanned by the knots
 */

 template <typename T, typename XkType, typename YkType, typename DydxkType,
          require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T, XkType, YkType, DydxkType> interpolate_1d_cubic_spline_eval(
    const T& x, const std::tuple<XkType, YkType, DydxkType>& W) {
  static constexpr const char* function = "interpolate_1d_cubic_spline_eval";
  const auto& xk = std::get<0>(W);
  const auto& yk = std::get<1>(W);
  const auto& dydxk = std::get<2>(W);

  check_size_match(function, "xk", xk.size(), "yk", yk.size());
  check_size_match(function, "xk", xk.size(), "dydxk", dydxk.size());
  check_finite(function, "x", x);

  const Eigen::Index K = xk.size();
  check_bounded(function, "x", x, xk.coeff(0), xk.coeff(K - 1));

  // Locate the interval i such that xk(i) <= x <= xk(i + 1)
  const auto* begin = xk.data();
  const auto* end = begin + K;
  const auto it = std::upper_bound(begin, end, x);
  Eigen::Index i = std::distance(begin, it) - 1;
  if (i == K - 1) {
    --i;
  }

  // Compute scaled parameter t and interval width h
  const auto h = xk.coeff(i + 1) - xk.coeff(i);
  const auto t = (x - xk.coeff(i)) / h;
  const auto t_minus_1 = t - 1.0;

  // Evaluate using standard cubic Hermite basis functions
  const auto y0 = yk.coeff(i);
  const auto y1 = yk.coeff(i + 1);
  const auto s0 = dydxk.coeff(i);
  const auto s1 = dydxk.coeff(i + 1);  
  
  const auto h0 = t_minus_1 * t_minus_1 * (y0 * (1.0 + 2.0 * t) + h * s0 * t);
  const auto h1 = t * t * (y1 * (3.0 - 2.0 * t) + h * s1 * t_minus_1);

  return h0 + h1;
}

// to be deleted 
// hermite interpolation function - wrapper for generic spline evaluation function 
// template <typename T, typename XkType, typename YkType, typename DydxkType>
// inline auto interpolate_1d_cubic_hermite_eval(const T& x, const std::tuple<XkType, YkType, DydxkType>& W) {
//   return interpolate_1d_cubic_spline_eval(x, W);
// }

// // Forward to your actual implementation
// template <typename EigVecX, typename EigVecY, typename EigVecDyDx>
// inline auto interpolate_1d_cubic_hermite_setup(const EigVecX& xk, const EigVecY& yk,
//                                       const EigVecDyDx& dydxk, std::ostream* pstream__) {
//   return interpolate_1d_cubic_hermite_setup(xk, yk, dydxk);
// }
// 
// template <typename T, typename XkType, typename YkType, typename DydxkType>
// inline auto interpolate_1d_cubic_hermite_eval(const T& x, const std::tuple<XkType, YkType, DydxkType>& W,
//                                      std::ostream* pstream__) {
//   return interpolate_1d_cubic_hermite_eval(x, W);
// }

}  // namespace math
}  // namespace stan

#endif
