#ifndef STAN_MATH_REV_META_EIGEN_PRODUCT_EVALUATOR_HOOKS_HPP
#define STAN_MATH_REV_META_EIGEN_PRODUCT_EVALUATOR_HOOKS_HPP

/**
 * @file
 * Specializations of Eigen's evaluators that compute products of matrices
 * with autodiff scalars with Stan's `multiply`, `elt_multiply`,
 * `scale_multiply` and `multiply_add`, so that `A * B` records one node on
 * the reverse-mode stack instead of one per scalar operation.
 *
 * Each specialization must be more specialized than Eigen's under partial
 * ordering. A condition (`enable_if_t`) can only go in an argument that
 * Eigen's specialization leaves free; where Eigen's has a computed default
 * or its own condition, ours must repeat it exactly or the two are
 * ambiguous. When upgrading Eigen, compare the specializations in
 * Eigen/src/Core/ProductEvaluators.h and CoreEvaluators.h with these.
 *
 * This file is included by stan/math/rev/meta.hpp and
 * stan/math/rev/core/Eigen_NumTraits.hpp, so every translation unit that
 * uses `var` sees it. The functions the specializations call are found by
 * argument-dependent lookup when a product is evaluated, so a translation
 * unit that evaluates a `var` product must include stan/math/rev.hpp.
 */

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <type_traits>

namespace stan {
namespace math {
namespace internal {

/**
 * True if either operand of an Eigen product has `var` scalars and neither
 * has complex scalars, so the product is computed by Stan's functions.
 * Products with a complex operand are left to Eigen: Stan's `multiply` for
 * them evaluates an Eigen `lazyProduct`, which would come back here.
 *
 * @tparam Lhs left operand expression type
 * @tparam Rhs right operand expression type
 */
template <typename Lhs, typename Rhs>
inline constexpr bool is_stan_ad_product_v
    = is_var<return_type_t<Lhs, Rhs>>::value;

template <typename T>
struct is_scalar_times_matrix_impl : std::false_type {};
template <typename S1, typename S2, typename Plain, typename Mat>
struct is_scalar_times_matrix_impl<Eigen::CwiseBinaryOp<
    Eigen::internal::scalar_product_op<S1, S2>,
    const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<S1>, Plain>,
    const Mat>> : std::true_type {};

    /**
 * True for the expression Eigen builds for a scalar times a matrix, `s * X`:
 * a `CwiseBinaryOp` of `scalar_product_op` whose left operand is a
 * `CwiseNullaryOp<scalar_constant_op>` holding `s`
 * (`xpr.lhs().functor().m_other`) and whose right operand is `X` (`xpr.rhs()`).
 *
 * @tparam T an Eigen expression type, with any reference or `const`
 */
template <typename T>
struct is_scalar_times_matrix : is_scalar_times_matrix_impl<std::decay_t<T>> {};

template <typename T>
struct is_matrix_times_scalar_impl : std::false_type {};
template <typename S1, typename S2, typename Mat, typename Plain>
struct is_matrix_times_scalar_impl<
    Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<S1, S2>, const Mat,
                         const Eigen::CwiseNullaryOp<
                             Eigen::internal::scalar_constant_op<S2>, Plain>>>
    : std::true_type {};

/**
 * True for the expression Eigen builds for a matrix times a scalar, `X * s`:
 * a `CwiseBinaryOp` of `scalar_product_op` whose left operand is `X`
 * (`xpr.lhs()`) and whose right operand is a
 * `CwiseNullaryOp<scalar_constant_op>` holding `s`
 * (`xpr.rhs().functor().m_other`).
 *
 * @tparam T an Eigen expression type, with any reference or `const`
 */
template <typename T>
struct is_matrix_times_scalar : is_matrix_times_scalar_impl<std::decay_t<T>> {};

template <typename T>
struct is_eigen_product_impl : std::false_type {};
template <typename Lhs, typename Rhs, int Options>
struct is_eigen_product_impl<Eigen::Product<Lhs, Rhs, Options>>
    : std::true_type {};

/**
 * True for an Eigen matrix product expression.
 *
 * @tparam T an Eigen expression type, with any reference or `const`
 */
template <typename T>
struct is_eigen_product : is_eigen_product_impl<std::decay_t<T>> {};

template <typename T>
struct is_eigen_cwise_product_impl : std::false_type {};
template <typename LhsScalar, typename RhsScalar, typename Lhs, typename Rhs>
struct is_eigen_cwise_product_impl<Eigen::CwiseBinaryOp<
    Eigen::internal::scalar_product_op<LhsScalar, RhsScalar>, Lhs, Rhs>>
    : std::true_type {};

/**
 * True for an Eigen coefficient-wise product expression: `A.cwiseProduct(B)`,
 * `a * b` of arrays, and `s * X` / `X * s`, which Eigen builds as a
 * coefficient-wise product with a constant expression.
 *
 * @tparam T an Eigen expression type, with any reference or `const`
 */
template <typename T>
struct is_eigen_cwise_product : is_eigen_cwise_product_impl<std::decay_t<T>> {};

/*! \brief Require type satisfies @ref is_eigen_product */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_product_t = require_t<is_eigen_product<T>>;

/*! \brief Require type satisfies @ref is_eigen_cwise_product */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_cwise_product_t = require_t<is_eigen_cwise_product<T>>;

/**
 * Empty base of every Stan specialization of an Eigen evaluator or
 * assignment below, so tests can check which specialization Eigen picks for
 * an expression (`std::is_base_of<eigen_ad_hook, evaluator<Xpr>>`).
 */
struct eigen_ad_hook {};

/**
 * `var` if Eigen's assignment of the product of `Lhs` and `Rhs` is computed
 * by Stan's functions (both operands are dense), otherwise no type. Used as
 * both scalars of the assignment functor in the `Assignment`
 * specializations below.
 *
 * @tparam Lhs left operand expression type
 * @tparam Rhs right operand expression type
 */
template <typename Lhs, typename Rhs>
using var_if_dense_product_t = std::enable_if_t<
    std::is_same<typename Eigen::internal::evaluator_traits<Lhs>::Shape,
                 Eigen::DenseShape>::value
        && std::is_same<typename Eigen::internal::evaluator_traits<Rhs>::Shape,
                        Eigen::DenseShape>::value,
    var>;

/**
 * Return the product of the operands of an Eigen matrix product. When an
 * operand is Eigen's `s * A` or `A * s` (Eigen rewrites `s * (A * B)` as
 * `(s * A) * B`) this is `scale_multiply(s, A, B)`, otherwise
 * `multiply(lhs, rhs)`.
 *
 * @tparam Lhs left operand expression type
 * @tparam Rhs right operand expression type
 * @param lhs left operand
 * @param rhs right operand
 */
template <typename Lhs, typename Rhs>
inline auto multiply_operands(Lhs&& lhs, Rhs&& rhs) {
  if constexpr (is_scalar_times_matrix<Lhs>::value) {
    // (s * A) * B
    return scale_multiply(std::forward<Lhs>(lhs).lhs().functor().m_other, std::forward<Lhs>(lhs).rhs(),
                          std::forward<Rhs>(rhs));
  } else if constexpr (is_matrix_times_scalar<Lhs>::value) {
    // (A * s) * B
    return scale_multiply(std::forward<Lhs>(lhs).rhs().functor().m_other, std::forward<Lhs>(lhs).lhs(),
                          std::forward<Rhs>(rhs));
  } else if constexpr (is_scalar_times_matrix<Rhs>::value) {
    // A * (s * B)
    return scale_multiply(std::forward<Rhs>(rhs).lhs().functor().m_other, std::forward<Lhs>(lhs),
                          std::forward<Rhs>(rhs).rhs());
  } else if constexpr (is_matrix_times_scalar<Rhs>::value) {
    // A * (B * s)
    return scale_multiply(std::forward<Rhs>(rhs).rhs().functor().m_other, std::forward<Lhs>(lhs),
                          std::forward<Rhs>(rhs).lhs());
  } else {
    return multiply(std::forward<Lhs>(lhs), std::forward<Rhs>(rhs));
  }
}

/**
 * Return the matrix product `xpr` computed by `multiply_operands`.
 *
 * @param xpr an Eigen matrix product
 */
template <typename XprType, require_eigen_product_t<XprType>* = nullptr>
inline auto evaluate_product(XprType&& xpr) {
  return multiply_operands(std::forward<XprType>(xpr).lhs(), std::forward<XprType>(xpr).rhs());
}

/**
 * Return the coefficient-wise product `xpr`: `scale_multiply(s, A, B)` for
 * `s * (A * B)` or `(A * B) * s`, `multiply(s, X)` for `s * X` or `X * s`,
 * and `elt_multiply(lhs, rhs)` otherwise. Stan's functions take matrices, so
 * array operands are viewed as matrices (`.matrix().derived()` is the
 * operand itself for a matrix).
 *
 * @param xpr an Eigen coefficient-wise product
 */
template <typename XprType, require_eigen_cwise_product_t<XprType>* = nullptr>
inline auto evaluate_product(XprType&& xpr) {
  // the operand types, without a reference or `const`
  using Xpr = std::decay_t<XprType>;
  if constexpr (is_scalar_times_matrix<XprType>::value) {
    if constexpr (is_eigen_product<typename Xpr::Rhs>::value) {
      // s * (A * B)
      return scale_multiply(std::forward<XprType>(xpr).lhs().functor().m_other, std::forward<XprType>(xpr).rhs().lhs(),
                            std::forward<XprType>(xpr).rhs().rhs());
    } else {
      // s * X
      return multiply(std::forward<XprType>(xpr).lhs().functor().m_other,
                      std::forward<XprType>(xpr).rhs().matrix().derived());
    }
  } else if constexpr (is_matrix_times_scalar<XprType>::value) {
    if constexpr (is_eigen_product<typename Xpr::Lhs>::value) {
      // (A * B) * s
      return scale_multiply(std::forward<XprType>(xpr).rhs().functor().m_other,
                             std::forward<XprType>(xpr).lhs().lhs(),
                             std::forward<XprType>(xpr).lhs().rhs());
    } else {
      // X * s
      return multiply(std::forward<XprType>(xpr).rhs().functor().m_other,
                      std::forward<XprType>(xpr).lhs().matrix().derived());
    }
  } else {
    return elt_multiply(std::forward<XprType>(xpr).lhs().matrix().derived(),
                        std::forward<XprType>(xpr).rhs().matrix().derived());
  }
}

/**
 * Body of the `product_evaluator` and `binary_evaluator` specializations
 * below. Like Eigen's product evaluator it holds the product, computed by
 * `evaluate_product`, and evaluates to it.
 *
 * @tparam XprType an Eigen matrix or coefficient-wise product
 */
template <typename XprType>
struct ad_product_evaluator
    : public Eigen::internal::evaluator<typename XprType::PlainObject>,
      eigen_ad_hook {
  using PlainObject = typename XprType::PlainObject;
  using Base = Eigen::internal::evaluator<PlainObject>;
  enum { Flags = Base::Flags | Eigen::EvalBeforeNestingBit };
  template <typename XprTypeT, stan::require_t<std::is_same<std::decay_t<XprTypeT>, std::decay_t<XprType>>>* = nullptr>
  explicit ad_product_evaluator(XprTypeT&& xpr) {
    if constexpr (is_stan_scalar<decltype(evaluate_product(xpr))>::value) {
      // a row vector times a column vector
      m_result.resize(1, 1);
      m_result.coeffRef(0, 0) = evaluate_product(std::forward<XprTypeT>(xpr));
    } else if constexpr (std::is_base_of<Eigen::ArrayBase<PlainObject>,
                                         PlainObject>::value) {
      m_result = evaluate_product(std::forward<XprTypeT>(xpr)).array();
    } else {
      m_result = evaluate_product(std::forward<XprTypeT>(xpr));
    }
    // Base points into m_result's storage, so it is constructed only once
    // m_result is filled: assigning a returned matrix can move a new buffer
    // into m_result.
    Eigen::internal::construct_at<Base>(this, m_result);
  }

 protected:
  PlainObject m_result;
};

}  // namespace internal
}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

/**
 * Dense-by-dense products where an operand has `var` scalars, computed by
 * `ad_product_evaluator` with `multiply_operands` (`multiply` or
 * `scale_multiply`).
 *
 * The `var` condition sits in the left shape argument: Eigen's generic
 * `product_evaluator` takes any type there, so this specialization is more
 * specialized than it and does not exist for other scalars. (Putting the
 * condition in the scalar arguments or in the type of `Options` makes the
 * two specializations ambiguous, because Eigen's generic one computes the
 * scalars and takes an `int` there.)
 *
 * Lazy products (`lazyProduct()`, and the `LazyCoeffBasedProductMode`
 * evaluator that derives from Eigen's lazy one) are left to Eigen: its lazy
 * evaluators fix both shapes, so this specialization would be ambiguous with
 * them.
 */
template <typename Lhs, typename Rhs, int Options, int ProductTag>
struct product_evaluator<
    Product<Lhs, Rhs, Options>, ProductTag,
    std::enable_if_t<stan::math::internal::is_stan_ad_product_v<
                         Lhs, Rhs> && Options != LazyProduct
                         && ProductTag != LazyCoeffBasedProductMode,
                     DenseShape>,
    DenseShape>
    : stan::math::internal::ad_product_evaluator<Product<Lhs, Rhs, Options>> {
  using XprType = Product<Lhs, Rhs, Options>;
  template <typename XprTypeT,
            stan::math::internal::require_eigen_product_t<XprTypeT>* = nullptr>
  explicit product_evaluator(XprTypeT&& xpr)
      : stan::math::internal::ad_product_evaluator<XprType>(
          std::forward<XprTypeT>(xpr)) {}
};

/**
 * Coefficient-wise products `lhs .* rhs` where an operand has `var` scalars
 * and neither is complex, computed by `ad_product_evaluator`.
 *
 * The condition sits in the functor argument of `CwiseBinaryOp`, which
 * Eigen's specialization leaves free, so this one is more specialized. The
 * functor is written as Eigen builds it, `scalar_product_op` of the two
 * operands' scalars, so `Lhs` and `Rhs` are deduced from the operands.
 */
template <typename Lhs, typename Rhs>
struct binary_evaluator<
    CwiseBinaryOp<
        std::enable_if_t<stan::math::internal::is_stan_ad_product_v<Lhs, Rhs>,
                         scalar_product_op<typename traits<Lhs>::Scalar,
                                           typename traits<Rhs>::Scalar>>,
        Lhs, Rhs>,
    IndexBased, IndexBased>
    : stan::math::internal::ad_product_evaluator<
          CwiseBinaryOp<scalar_product_op<typename traits<Lhs>::Scalar,
                                          typename traits<Rhs>::Scalar>,
                        Lhs, Rhs>> {
  using XprType = CwiseBinaryOp<scalar_product_op<typename traits<Lhs>::Scalar,
                                                  typename traits<Rhs>::Scalar>,
                                Lhs, Rhs>;
  explicit binary_evaluator(const XprType& xpr)
      : stan::math::internal::ad_product_evaluator<XprType>(xpr) {}
};

/**
 * Eigen assigns a product to a destination (`dst = A * B`, and
 * `dst.noalias() += A * B` / `-=`) with its own `Assignment`
 * specializations, which compute the product without `product_evaluator`.
 * Eigen has one specialization per assignment functor, so these are three
 * as well. The condition (`var` product of dense operands) sits in both
 * scalars of the functor, which Eigen's specializations leave free, so these
 * are more specialized; other products run Eigen's own specialization. The
 * last argument must repeat Eigen's own condition token for token (a
 * different condition makes the two ambiguous).
 */
template <typename DstXprType, typename Lhs, typename Rhs, int Options>
struct Assignment<
    DstXprType, Product<Lhs, Rhs, Options>,
    assign_op<stan::math::internal::var_if_dense_product_t<Lhs, Rhs>,
              stan::math::internal::var_if_dense_product_t<Lhs, Rhs>>,
    Dense2Dense,
    std::enable_if_t<(Options == DefaultProduct
                      || Options == AliasFreeProduct)>>
    : stan::math::internal::eigen_ad_hook {
  template <typename Func, typename SrcXpr,
            stan::math::internal::require_eigen_product_t<SrcXpr>* = nullptr>
  static void run(DstXprType& dst, SrcXpr&& src, Func&& /* func */) {
    using stan::math::internal::multiply_operands;
    if constexpr (stan::is_stan_scalar<decltype(
                      multiply_operands(src.lhs(), src.rhs()))>::value) {
      // a row vector times a column vector
      if (dst.rows() != 1 || dst.cols() != 1) {
        dst.resize(1, 1);
      }
      dst.coeffRef(0, 0) = multiply_operands(std::forward<SrcXpr>(src).lhs(), std::forward<SrcXpr>(src).rhs());
    } else {
      dst = multiply_operands(std::forward<SrcXpr>(src).lhs(), std::forward<SrcXpr>(src).rhs());
    }
  }
};

template <typename DstXprType, typename Lhs, typename Rhs, int Options>
struct Assignment<
    DstXprType, Product<Lhs, Rhs, Options>,
    add_assign_op<stan::math::internal::var_if_dense_product_t<Lhs, Rhs>,
                  stan::math::internal::var_if_dense_product_t<Lhs, Rhs>>,
    Dense2Dense,
    std::enable_if_t<(Options == DefaultProduct
                      || Options == AliasFreeProduct)>>
    : stan::math::internal::eigen_ad_hook {
  template <typename Func, typename SrcType, stan::math::internal::require_eigen_product_t<SrcType>* = nullptr>
  static void run(DstXprType& dst, SrcType&& src, Func&& /* func */) {
    dst = multiply_add(std::forward<SrcType>(src).lhs(), std::forward<SrcType>(src).rhs(), dst);
  }
};

template <typename DstXprType, typename Lhs, typename Rhs, int Options>
struct Assignment<
    DstXprType, Product<Lhs, Rhs, Options>,
    sub_assign_op<stan::math::internal::var_if_dense_product_t<Lhs, Rhs>,
                  stan::math::internal::var_if_dense_product_t<Lhs, Rhs>>,
    Dense2Dense,
    std::enable_if_t<(Options == DefaultProduct
                      || Options == AliasFreeProduct)>>
    : stan::math::internal::eigen_ad_hook {
  template <typename Func, typename SrcType, stan::math::internal::require_eigen_product_t<SrcType>* = nullptr>
  static void run(DstXprType& dst, SrcType&& src, Func&& /* func */) {
    dst = subtract(
        dst, stan::math::internal::multiply_operands(std::forward<SrcType>(src).lhs(), std::forward<SrcType>(src).rhs()));
  }
};

}  // namespace internal
}  // namespace Eigen

#endif
