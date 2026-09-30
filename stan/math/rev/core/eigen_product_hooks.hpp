#ifndef STAN_MATH_REV_CORE_EIGEN_PRODUCT_HOOKS_HPP
#define STAN_MATH_REV_CORE_EIGEN_PRODUCT_HOOKS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/eigen_product_kernel.hpp>
#include <stan/math/rev/core/std_complex.hpp>
#include <stan/math/rev/core/var.hpp>
#include <complex>
#include <type_traits>

/**
 * @file
 * Routes Eigen products whose operands have `var` or `std::complex<var>`
 * scalars (mixed with arithmetic scalars) to the reverse-mode kernels in
 * stan/math/rev/core/eigen_product_kernel.hpp.
 *
 *   - matrix products `A * B` (every Eigen product tag, `lazyProduct`,
 *     `noalias()`, `=`, `+=` and `-=`) through the `product_evaluator` and
 *     `Assignment` hooks of
 * stan/math/prim/meta/eigen_product_evaluator_hooks.hpp;
 *   - coefficient-wise products `a.array() * b.array()`, `cwiseProduct`,
 *     and products with a scalar `s * A`, `A * s` through `binary_evaluator`.
 *
 * Scalar factors are folded into the product they scale, so `s * (A * B)`,
 * `(A * B) * s` and `(s * A) * B` each put one callback on the stack.
 *
 * The product is evaluated in full when its evaluator is created, even when
 * only some coefficients are read afterwards (e.g. `(a * b).head(2)`). An
 * unread coefficient's adjoint is zero, but the callback still multiplies it
 * by the operand values, so a non-finite operand value can turn it into NaN.
 * In particular, a product inside a branch of `.select()` is evaluated for
 * every coefficient, unlike other expressions there.
 */

namespace stan {
namespace math {
namespace internal {

template <typename Lhs, typename Rhs>
inline constexpr bool is_rev_product_v
    = is_rev_product_operand_scalar<
          typename Eigen::internal::traits<Lhs>::Scalar>::value&&
          is_rev_product_operand_scalar<
              typename Eigen::internal::traits<Rhs>::Scalar>::value
      && (is_rev_product_scalar<
              typename Eigen::internal::traits<Lhs>::Scalar>::value
          || is_rev_product_scalar<
              typename Eigen::internal::traits<Rhs>::Scalar>::value);

/**
 * Reverse-mode matrix product.
 */
template <typename Lhs, typename Rhs>
struct eigen_matrix_product_impl<Lhs, Rhs,
                                 std::enable_if_t<is_rev_product_v<Lhs, Rhs>>> {
  using result_t = typename Eigen::Product<Lhs, Rhs>::PlainObject;
  static result_t run(const Lhs& lhs, const Rhs& rhs) {
    return rev_matrix_product<result_t>(lhs, rhs);
  }
};

/**
 * Reverse-mode coefficient-wise product, or product with a scalar.
 */
template <typename S1, typename S2, typename Lhs, typename Rhs>
struct eigen_cwise_product_impl<
    Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<S1, S2>, Lhs, Rhs>,
    std::enable_if_t<is_rev_product_v<Lhs, Rhs>>> {
  using xpr_t = Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<S1, S2>,
                                     Lhs, Rhs>;
  using result_t = typename xpr_t::PlainObject;
  static result_t run(const xpr_t& xpr) {
    return rev_cwise_product<result_t>(xpr);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

#define STAN_REV_PRODUCT_NO_TPARAMS()
#define STAN_REV_PRODUCT_S_TPARAMS() typename S_,
STAN_EIGEN_AD_SCALAR_PATTERNS(STAN_EIGEN_AD_CWISE_PRODUCT_EVALUATOR,
                              STAN_REV_PRODUCT_NO_TPARAMS,
                              STAN_REV_PRODUCT_S_TPARAMS, stan::math::var)
STAN_EIGEN_AD_PRODUCT_ASSIGNMENT(STAN_REV_PRODUCT_NO_TPARAMS, stan::math::var)
STAN_EIGEN_AD_PRODUCT_ASSIGNMENT(STAN_REV_PRODUCT_NO_TPARAMS,
                                 std::complex<stan::math::var>)
#undef STAN_REV_PRODUCT_NO_TPARAMS
#undef STAN_REV_PRODUCT_S_TPARAMS

}  // namespace internal
}  // namespace Eigen

#endif
