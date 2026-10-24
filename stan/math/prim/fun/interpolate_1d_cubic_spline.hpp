#ifndef STAN_MATH_PRIM_FUN_INTERPOLATE_1D_CUBIC_SPLINE_HPP
#define STAN_MATH_PRIM_FUN_INTERPOLATE_1D_CUBIC_SPLINE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/meta.hpp>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/interpolators/makima.hpp>

namespace stan {
namespace math {

/**
 * Storage tuple type for 1D cubic spline interpolants.
 *
 * The tuple holds the knot locations, knot values, and derivative values as
 * owning Eigen column vectors, used as the precomputed representation consumed
 * by `interpolate_1d_cubic_spline`.
 *
 * @tparam XkScalar scalar type of the knot locations
 * @tparam YkScalar scalar type of the knot values
 * @tparam DydxkScalar scalar type of the derivative values
 */
template <typename XkScalar, typename YkScalar, typename DydxkScalar>
using interpolate_1d_cubic_spline_store_t
    = std::tuple<Eigen::Matrix<XkScalar, Eigen::Dynamic, 1>,
                 Eigen::Matrix<YkScalar, Eigen::Dynamic, 1>,
                 Eigen::Matrix<DydxkScalar, Eigen::Dynamic, 1>>;

/**
 * Precompute storage for cubic Hermite interpolation.
 *
 * This function validates the input knot locations, values, and derivatives,
 * then stores them in a tuple for use with `interpolate_1d_cubic_spline`.
 *
 * @note Because the derivatives at the knots are supplied directly (rather
 *   than derived by branching on `yk` as in PCHIP/MAKIMA), the interpolated
 *   curve is smooth and the resulting log density is differentiable
 *   everywhere with respect to `yk`.
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
inline interpolate_1d_cubic_spline_store_t<return_type_t<EigVecX>,
                                           return_type_t<EigVecY>,
                                           return_type_t<EigVecDyDx>>
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
 * @param W Tuple of (xk, yk, dydxk) from one of the `_setup` functions.
 * @return The interpolated value at `x`.
 * @throw std::invalid_argument if `W` does not have the expected structure
 * @throw std::domain_error if `x` is not finite or is outside of the range
 *   spanned by the knots
 */
 template <typename T, typename XkType, typename YkType, typename DydxkType,
          require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T, XkType, YkType, DydxkType> interpolate_1d_cubic_spline(
    const T& x, const std::tuple<XkType, YkType, DydxkType>& W) {
  static constexpr const char* function = "interpolate_1d_cubic_spline";
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

/**
 * Precompute storage for PCHIP interpolation.
 *
 * Computes PCHIP (monotone-preserving) derivatives at the knots via Boost,
 * then delegates to `interpolate_1d_cubic_hermite_setup` so the result can be
 * evaluated with `interpolate_1d_cubic_spline`. The interpolation method
 * itself -- including the formulas for the end slopes and the requirement of
 * at least four points -- is defined by Boost; see the documentation for
 * `boost::math::interpolators::pchip`.
 *
 * @note The PCHIP end slopes computed by Boost differ from those used by scipy
 *   and MATLAB, so results at and near the boundary knots may not match those
 *   libraries.
 *
 * @note The PCHIP slope formulas branch on the data (in `yk`), so the
 *   interpolated curve has kinks in `yk` and the resulting log density is not
 *   differentiable everywhere with respect to `yk`. If differentiability with
 *   respect to `yk` is required, use `interpolate_1d_cubic_hermite_setup`
 *   instead, which takes known derivatives and has no such problem.
 *
 * @tparam EigVecX type of the knot locations
 * @tparam EigVecY type of the knot values
 * @param xk Strictly increasing knot locations (at least four).
 * @param yk Interpolation values at the knots.
 * @return A `(xk, yk, dydxk)` tuple with PCHIP derivatives.
 * @throw std::invalid_argument if `xk` and `yk` have different sizes
 * @throw std::domain_error if there are fewer than four knots, any input is
 *   not finite, or `xk` is not strictly increasing
 */
template <typename EigVecX, typename EigVecY,
          require_all_eigen_col_vector_t<EigVecX, EigVecY>* = nullptr>
inline interpolate_1d_cubic_spline_store_t<return_type_t<EigVecX>,
                                           return_type_t<EigVecY>,
                                           return_type_t<EigVecX, EigVecY>>
interpolate_1d_pchip_setup(const EigVecX& xk, const EigVecY& yk) {
  static constexpr const char* function = "interpolate_1d_pchip_setup";
  const auto& xk_ref = to_ref(xk);
  const auto& yk_ref = to_ref(yk);

  check_greater_or_equal(function, "number of knots", xk_ref.size(), 4);
  check_size_match(function, "xk", xk_ref.size(), "yk", yk_ref.size());
  check_finite(function, "xk", xk_ref);
  check_finite(function, "yk", yk_ref);
  check_ordered(function, "xk", xk_ref);

  using scalar_t = return_type_t<EigVecX, EigVecY>;

  // Boost takes std::vector containers
  std::vector<scalar_t> xk_vec(xk_ref.data(), xk_ref.data() + xk_ref.size());
  std::vector<scalar_t> yk_vec(yk_ref.data(), yk_ref.data() + yk_ref.size());
  boost::math::interpolators::pchip<std::vector<scalar_t>> pchip_interpolator(
      std::move(xk_vec), std::move(yk_vec));

  // PCHIP derivative at each knot
  std::vector<scalar_t> dydxk_vec;
  dydxk_vec.reserve(xk_ref.size());
  for (std::size_t i = 0; i < xk_ref.size(); ++i) {
    dydxk_vec.push_back(pchip_interpolator.prime(xk_ref[i]));
  }
  Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> dydxk_out
      = Eigen::Map<Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>>(
          dydxk_vec.data(), dydxk_vec.size());

  // Reuse the cubic Hermite storage/eval path
  return interpolate_1d_cubic_hermite_setup(xk_ref, yk_ref, dydxk_out);
}

/**
 * Precompute storage for MAKIMA interpolation.
 *
 * Computes modified Akima derivatives at the knots via Boost, then delegates to
 * `interpolate_1d_cubic_hermite_setup` so the result can be evaluated with
 * `interpolate_1d_cubic_spline`. The interpolation method itself -- including
 * the formulas for the end slopes and the requirement of at least four points
 * -- is defined by Boost; see the documentation for
 * `boost::math::interpolators::makima`.
 *
 * @note The MAKIMA slope formulas branch on the data (in `yk`), so the
 *   interpolated curve has kinks in `yk` and the resulting log density is not
 *   differentiable everywhere with respect to `yk`. If differentiability with
 *   respect to `yk` is required, use `interpolate_1d_cubic_hermite_setup`
 *   instead, which takes known derivatives and has no such problem.
 *
 * @tparam EigVecX type of the knot locations
 * @tparam EigVecY type of the knot values
 * @param xk Strictly increasing knot locations (at least four).
 * @param yk Interpolation values at the knots.
 * @return A `(xk, yk, dydxk)` tuple with MAKIMA derivatives.
 * @throw std::invalid_argument if `xk` and `yk` have different sizes
 * @throw std::domain_error if there are fewer than four knots, any input is
 *   not finite, or `xk` is not strictly increasing
 */
template <typename EigVecX, typename EigVecY,
          require_all_eigen_col_vector_t<EigVecX, EigVecY>* = nullptr>
inline interpolate_1d_cubic_spline_store_t<return_type_t<EigVecX>,
                                           return_type_t<EigVecY>,
                                           return_type_t<EigVecX, EigVecY>>
interpolate_1d_makima_setup(const EigVecX& xk, const EigVecY& yk) {
  static constexpr const char* function = "interpolate_1d_makima_setup";
  const auto& xk_ref = to_ref(xk);
  const auto& yk_ref = to_ref(yk);

  check_greater_or_equal(function, "number of knots", xk_ref.size(), 4);
  check_size_match(function, "xk", xk_ref.size(), "yk", yk_ref.size());
  check_finite(function, "xk", xk_ref);
  check_finite(function, "yk", yk_ref);
  check_ordered(function, "xk", xk_ref);

  using scalar_t = return_type_t<EigVecX, EigVecY>;

  // Boost takes std::vector containers
  std::vector<scalar_t> xk_vec(xk_ref.data(), xk_ref.data() + xk_ref.size());
  std::vector<scalar_t> yk_vec(yk_ref.data(), yk_ref.data() + yk_ref.size());
  boost::math::interpolators::makima<std::vector<scalar_t>> makima_interpolator(
      std::move(xk_vec), std::move(yk_vec));

  // MAKIMA derivative at each knot
  std::vector<scalar_t> dydxk_vec;
  dydxk_vec.reserve(xk_ref.size());
  for (std::size_t i = 0; i < xk_ref.size(); ++i) {
    dydxk_vec.push_back(makima_interpolator.prime(xk_ref[i]));
  }
  Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> dydxk_out
      = Eigen::Map<Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>>(
          dydxk_vec.data(), dydxk_vec.size());

  // Reuse the cubic Hermite storage/eval path
  return interpolate_1d_cubic_hermite_setup(xk_ref, yk_ref, dydxk_out);
}

}  // namespace math
}  // namespace stan

#endif

