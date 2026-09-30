#ifndef STAN_MATH_FWD_CORE_EIGEN_PRODUCT_HOOKS_HPP
#define STAN_MATH_FWD_CORE_EIGEN_PRODUCT_HOOKS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core/eigen_product_kernel.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/core/std_complex.hpp>
#include <complex>
#include <type_traits>

/**
 * @file
 * Routes Eigen products whose operands have `fvar<T>` or
 * `std::complex<fvar<T>>` scalars (mixed with arithmetic scalars) to the
 * forward-mode kernels in stan/math/fwd/core/eigen_product_kernel.hpp, the
 * same way stan/math/rev/core/eigen_product_hooks.hpp does for `var`.
 */

namespace stan {
namespace math {
namespace internal {

template <typename S, typename T>
struct is_fwd_product_operand_scalar
    : std::integral_constant<bool, std::is_arithmetic<S>::value> {};
template <typename S, typename T>
struct is_fwd_product_operand_scalar<std::complex<S>, T>
    : std::integral_constant<bool, std::is_arithmetic<S>::value> {};
template <typename T>
struct is_fwd_product_operand_scalar<fvar<T>, T> : std::true_type {};
template <typename T>
struct is_fwd_product_operand_scalar<std::complex<fvar<T>>, T>
    : std::true_type {};

template <typename LS, typename RS, typename T>
inline constexpr bool is_fwd_product_scalars_v
    = !std::is_void<T>::value && is_fwd_product_operand_scalar<LS, T>::value
      && is_fwd_product_operand_scalar<RS, T>::value;

/**
 * The inner scalar `T` shared by the `fvar` operands of a product, or `void`
 * if neither operand is an `fvar`.
 */
template <typename LS, typename RS>
using fwd_product_inner_of_t
    = std::conditional_t<std::is_void<fwd_product_inner_t<LS>>::value,
                         fwd_product_inner_t<RS>, fwd_product_inner_t<LS>>;

template <typename Lhs, typename Rhs>
inline constexpr bool is_fwd_product_v = is_fwd_product_scalars_v<
    typename Eigen::internal::traits<Lhs>::Scalar,
    typename Eigen::internal::traits<Rhs>::Scalar,
    fwd_product_inner_of_t<typename Eigen::internal::traits<Lhs>::Scalar,
                           typename Eigen::internal::traits<Rhs>::Scalar>>;

/**
 * Forward-mode matrix product.
 */
template <typename Lhs, typename Rhs>
struct eigen_matrix_product_impl<Lhs, Rhs,
                                 std::enable_if_t<is_fwd_product_v<Lhs, Rhs>>> {
  using result_t = typename Eigen::Product<Lhs, Rhs>::PlainObject;
  static result_t run(const Lhs& lhs, const Rhs& rhs) {
    return fwd_matrix_product<result_t>(lhs, rhs);
  }
};

/**
 * Forward-mode coefficient-wise product, or product with a scalar.
 */
template <typename S1, typename S2, typename Lhs, typename Rhs>
struct eigen_cwise_product_impl<
    Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<S1, S2>, Lhs, Rhs>,
    std::enable_if_t<is_fwd_product_v<Lhs, Rhs>>> {
  using xpr_t = Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<S1, S2>,
                                     Lhs, Rhs>;
  using result_t = typename xpr_t::PlainObject;
  static auto run(const xpr_t& xpr) { return fwd_cwise_product<result_t>(xpr); }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

#define STAN_FWD_PRODUCT_T_TPARAMS() typename T_,
#define STAN_FWD_PRODUCT_TS_TPARAMS() typename T_, typename S_,
STAN_EIGEN_AD_SCALAR_PATTERNS(STAN_EIGEN_AD_CWISE_PRODUCT_EVALUATOR,
                              STAN_FWD_PRODUCT_T_TPARAMS,
                              STAN_FWD_PRODUCT_TS_TPARAMS, stan::math::fvar<T_>)
STAN_EIGEN_AD_PRODUCT_ASSIGNMENT(STAN_FWD_PRODUCT_T_TPARAMS,
                                 stan::math::fvar<T_>)
STAN_EIGEN_AD_PRODUCT_ASSIGNMENT(STAN_FWD_PRODUCT_T_TPARAMS,
                                 std::complex<stan::math::fvar<T_>>)
#undef STAN_FWD_PRODUCT_T_TPARAMS
#undef STAN_FWD_PRODUCT_TS_TPARAMS

}  // namespace internal
}  // namespace Eigen

#endif
