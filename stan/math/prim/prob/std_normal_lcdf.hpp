#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/std_normal_lcdf_impl.hpp>

namespace stan {
namespace math {
namespace internal {
constexpr char std_normal_lcdf_func[] = "std_normal_lcdf";
}  // namespace internal

/** \ingroup prob_dists
 * @brief Calculates the log of the cdf of the standard normal distribution
 *
 * The piecewise value and derivative live in `internal::std_normal_lcdf_value`
 * and `internal::std_normal_lcdf_grad` in prim/prob/std_normal_lcdf_impl.hpp,
 * so that other distributions can share them. Their structure is the same as
 * `normal_lcdf`, and the cutoffs, their provenance and the measurements
 * behind them are documented once, on that function. See
 * prim/prob/normal_lcdf.hpp. That covers the A&S 7.1.26, Cody (1969) and
 * DLMF 7.12.1 references, the R `pnorm` and SciPy `log_ndtr`
 * cross-references, why the erfc/Cody crossover sits at 4, the
 * stan-dev/math#1411 origin of the interior Taylor cutoffs, and the two
 * cutoff tables.
 *
 * Two differences apply when reading it here. The scaled variable is
 * `scaled_y = y * INV_SQRT_TWO`, not `scaled_diff = (y - mu) / (sigma *
 * SQRT_TWO)`; since `mu = 0` and `sigma = 1` the two coincide, so every
 * cutoff value transfers unchanged. And the test that enforces the
 * `worst in-range` column for this function is the branch_accuracy test in
 * mix/prob/std_normal_cdf_log_test.cpp.
 *
 * @tparam func name reported by the error checks. Reflected distributions
 *   such as `std_normal_lccdf` delegate here and pass their own name so that
 *   exceptions name the function the user actually called.
 * @tparam T_y A vector or scalar type for the random variable.
 * @param y (Sequence of) scalar(s).
 * @return The log of the standard normal cdf evaluated at the specified
 *   argument. If given a container, the log of the product of the cdfs.
 */
template <
    const char* func = internal::std_normal_lcdf_func, typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_lcdf(const T_y& y) {
  using T_partials_return = partials_return_t<T_y>;
  using T_y_ref = ref_type_t<T_y>;
  static constexpr const char* function = func;
  T_y_ref y_ref = y;
  check_not_nan(function, "Random variable", y_ref);

  if (size_zero(y)) {
    return 0;
  }

  T_partials_return lcdf(0.0);
  auto ops_partials = make_partials_propagator(y_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  size_t N = stan::math::size(y);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return scaled_y = y_dbl * INV_SQRT_TWO;

    lcdf += internal::std_normal_lcdf_value(scaled_y);
    // The upper branch is the only one that can carry a non-finite value
    // through from a finite input, so it is guarded exactly where it was.
    if (scaled_y > 0.0 && !is_not_nan(lcdf)) {
      lcdf = 0;
    }

    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials)[n]
          += internal::std_normal_lcdf_grad(scaled_y) * INV_SQRT_TWO;
    }
  }

  return ops_partials.build(lcdf);
}

}  // namespace math
}  // namespace stan
#endif
