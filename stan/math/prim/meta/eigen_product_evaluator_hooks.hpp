#ifndef STAN_MATH_PRIM_META_EIGEN_PRODUCT_EVALUATOR_HOOKS_HPP
#define STAN_MATH_PRIM_META_EIGEN_PRODUCT_EVALUATOR_HOOKS_HPP

/**
 * @file
 * Hooks that let Stan evaluate Eigen matrix products whose operands have
 * autodiff scalars.
 *
 * Eigen dispatches `lhs * rhs` for dense operands through
 * `Eigen::internal::product_evaluator`. Its scalar template parameters are
 * defaulted from the operand types, so a partial specialization that fixes
 * them is ambiguous with Eigen's generic specialization. Instead, the
 * specializations below match every dense-by-dense product and choose, from
 * the operand scalar types, between
 *
 *   - `stan::math::internal::ad_product_evaluator`, which forwards to
 *     `eigen_matrix_product_impl<Lhs, Rhs>::run(lhs, rhs)`, a customization
 *     point defined by the autodiff layers (`rev`, `fwd`); and
 *   - Eigen's own evaluation for all other scalars (e.g. `double`), which is
 *     reproduced exactly: `eigen_dense_product_evaluator` is a copy of the
 *     generic `product_evaluator` in Eigen/src/Core/ProductEvaluators.h, and
 *     lazy products derive from Eigen's coefficient-based evaluator.
 *
 * Because these specializations are selected for every dense product, they
 * are included from stan/math/prim/fun/Eigen.hpp so that every translation
 * unit sees the same definition. When upgrading Eigen, check the generic
 * `product_evaluator` and the `LazyProduct` evaluator in
 * Eigen/src/Core/ProductEvaluators.h against the copies here.
 */

#include <stan/math/prim/fun/Eigen.hpp>
#include <complex>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * True for scalar types whose Eigen products Stan evaluates itself. The
 * primary template is false; `rev` specializes it for `var` and `fwd` for
 * `fvar<T>`. Complex scalars inherit from their value type.
 *
 * @tparam Scalar scalar type of an Eigen expression
 */
template <typename Scalar>
struct is_stan_ad_product_scalar : std::false_type {};

template <typename T>
struct is_stan_ad_product_scalar<std::complex<T>>
    : is_stan_ad_product_scalar<T> {};

/**
 * True if either operand of an Eigen product has a scalar type for which
 * Stan evaluates the product.
 *
 * @tparam Lhs left operand expression type
 * @tparam Rhs right operand expression type
 */
template <typename Lhs, typename Rhs>
inline constexpr bool is_stan_ad_product_v
    = is_stan_ad_product_scalar<
          typename Eigen::internal::traits<Lhs>::Scalar>::value
      || is_stan_ad_product_scalar<
          typename Eigen::internal::traits<Rhs>::Scalar>::value;

/**
 * Customization point for dense matrix products with autodiff operands.
 * Specializations provide `static auto run(const Lhs&, const Rhs&)`, which
 * returns an object with an Eigen evaluator holding the product, usually the
 * product's plain object type. The primary template is intentionally left
 * undefined so that a product with an autodiff scalar and no implementation
 * is a compile error.
 *
 * @tparam Lhs left operand expression type
 * @tparam Rhs right operand expression type
 * @tparam Enable used by specializations for SFINAE
 */
template <typename Lhs, typename Rhs, typename Enable = void>
struct eigen_matrix_product_impl;

/**
 * Customization point for coefficient-wise products (including the product
 * of a matrix and a scalar) with autodiff operands. Specializations provide
 * `static auto run(const Xpr&)` for the `CwiseBinaryOp` expression `Xpr`,
 * returning either an evaluated plain object or a lazy expression; the
 * implementation decides which.
 *
 * @tparam Xpr the `CwiseBinaryOp` of a `scalar_product_op`
 * @tparam Enable used by specializations for SFINAE
 */
template <typename Xpr, typename Enable = void>
struct eigen_cwise_product_impl;

/**
 * Owns the object returned by a product customization point. It is a base
 * class of `forwarding_evaluator` listed before the evaluator base, so the
 * result is constructed before the evaluator that refers to it.
 *
 * @tparam Result type returned by the customization point
 */
template <typename Result>
struct forwarded_result {
  Result result_;
  template <typename F>
  explicit forwarded_result(F&& make) : result_(make()) {}
};

/**
 * Eigen evaluator for an expression whose value is computed by a Stan
 * customization point. The returned object is stored and read through its
 * own Eigen evaluator. When the returned object is a plain matrix or array
 * the evaluator sets `EvalBeforeNestingBit`, so an enclosing expression
 * evaluates this one once instead of once per coefficient access.
 *
 * Eigen never copies evaluators (`evaluator_base` is noncopyable), so the
 * inner evaluator's pointer into `result_` stays valid.
 *
 * @tparam Result type returned by the customization point
 */
template <typename Result>
struct forwarding_evaluator : private forwarded_result<Result>,
                              public Eigen::internal::evaluator<Result> {
  using holder_t = forwarded_result<Result>;
  using Base = Eigen::internal::evaluator<Result>;
  static constexpr bool is_plain_result
      = std::is_base_of<Eigen::PlainObjectBase<Result>, Result>::value;
  enum {
    Flags
    = static_cast<int>(Base::Flags)
      | (is_plain_result ? static_cast<int>(Eigen::EvalBeforeNestingBit) : 0)
  };

  template <typename F>
  explicit forwarding_evaluator(F&& make)
      : holder_t(std::forward<F>(make)), Base(this->holder_t::result_) {}
};

template <typename Lhs, typename Rhs>
using matrix_product_result_t
    = std::decay_t<decltype(eigen_matrix_product_impl<Lhs, Rhs>::run(
        std::declval<const Lhs&>(), std::declval<const Rhs&>()))>;

/**
 * Evaluator of a dense product with an autodiff operand. It forwards to
 * `eigen_matrix_product_impl`.
 *
 * @tparam Lhs left operand expression type
 * @tparam Rhs right operand expression type
 * @tparam Options Eigen product options (default, lazy, alias free)
 */
template <typename Lhs, typename Rhs, int Options>
struct ad_product_evaluator
    : forwarding_evaluator<matrix_product_result_t<Lhs, Rhs>> {
  using XprType = Eigen::Product<Lhs, Rhs, Options>;
  using Base = forwarding_evaluator<matrix_product_result_t<Lhs, Rhs>>;
  static constexpr bool is_stan_ad_product_hook = true;
  explicit ad_product_evaluator(const XprType& xpr)
      : Base([&xpr]() {
          return eigen_matrix_product_impl<Lhs, Rhs>::run(xpr.lhs(), xpr.rhs());
        }) {}
};

/**
 * Evaluator of a coefficient-wise product with an autodiff operand. It
 * forwards to `eigen_cwise_product_impl`.
 *
 * @tparam Xpr the `CwiseBinaryOp` of a `scalar_product_op`
 */
template <typename Xpr>
struct ad_cwise_product_evaluator
    : forwarding_evaluator<std::decay_t<decltype(
          eigen_cwise_product_impl<Xpr>::run(std::declval<const Xpr&>()))>> {
  using XprType = Xpr;
  using Base = forwarding_evaluator<std::decay_t<decltype(
      eigen_cwise_product_impl<Xpr>::run(std::declval<const Xpr&>()))>>;
  static constexpr bool is_stan_ad_product_hook = true;
  explicit ad_cwise_product_evaluator(const Xpr& xpr)
      : Base([&xpr]() { return eigen_cwise_product_impl<Xpr>::run(xpr); }) {}
};

/**
 * Detects the scalar-times-expression shapes Eigen builds for `s * X` and
 * `X * s` (see `EIGEN_SCALAR_BINARYOP_EXPR_RETURN_TYPE` in Eigen's
 * util/Macros.h), so a kernel can fold the scalar into one callback.
 */
template <typename T>
struct scaled_expression : std::false_type {};

template <typename S1, typename S2, typename Plain, typename Inner>
struct scaled_expression<Eigen::CwiseBinaryOp<
    Eigen::internal::scalar_product_op<S1, S2>,
    const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<S1>, Plain>,
    const Inner>> : std::true_type {
  template <typename X>
  static const S1& scale(const X& x) {
    return x.lhs().functor().m_other;
  }
  template <typename X>
  static const Inner& inner(const X& x) {
    return x.rhs();
  }
};

template <typename S1, typename S2, typename Plain, typename Inner>
struct scaled_expression<Eigen::CwiseBinaryOp<
    Eigen::internal::scalar_product_op<S1, S2>, const Inner,
    const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<S2>,
                                Plain>>> : std::true_type {
  template <typename X>
  static const S2& scale(const X& x) {
    return x.rhs().functor().m_other;
  }
  template <typename X>
  static const Inner& inner(const X& x) {
    return x.lhs();
  }
};

template <typename S1, typename S2, typename Plain1, typename Plain2>
struct scaled_expression<
    Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<S1, S2>,
                         const Eigen::CwiseNullaryOp<
                             Eigen::internal::scalar_constant_op<S1>, Plain1>,
                         const Eigen::CwiseNullaryOp<
                             Eigen::internal::scalar_constant_op<S2>, Plain2>>>
    : std::true_type {
  using Inner
      = const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<S2>,
                                    Plain2>;
  template <typename X>
  static const S1& scale(const X& x) {
    return x.lhs().functor().m_other;
  }
  template <typename X>
  static Inner& inner(const X& x) {
    return x.rhs();
  }
};

/**
 * Maximum number of nested scalar factors folded out of one operand. Deeper
 * nests are left in the operand, which then goes through its own product
 * kernel.
 */
constexpr int max_peel_depth = 4;

/**
 * Strip up to `max_peel_depth` scalar factors off `x` into `scales` and
 * return the remaining expression.
 */
template <int Depth = 0, typename Scales, typename X>
inline decltype(auto) peel_scales(Scales& scales, const X& x) {
  if constexpr (scaled_expression<X>::value && Depth < max_peel_depth) {
    scales.push(scaled_expression<X>::scale(x));
    return peel_scales<Depth + 1>(scales, scaled_expression<X>::inner(x));
  } else {
    return static_cast<const X&>(x);
  }
}

/**
 * True for an Eigen product of two dense operands.
 */
template <typename T>
struct is_dense_product : std::false_type {};
template <typename Lhs, typename Rhs, int Options>
struct is_dense_product<Eigen::Product<Lhs, Rhs, Options>>
    : std::integral_constant<
          bool,
          std::is_same<typename Eigen::internal::evaluator_traits<Lhs>::Shape,
                       Eigen::DenseShape>::value
              && std::is_same<
                  typename Eigen::internal::evaluator_traits<Rhs>::Shape,
                  Eigen::DenseShape>::value> {};

/**
 * True if the Eigen evaluator `Evaluator` is one of Stan's autodiff product
 * evaluators, i.e. the expression it evaluates is computed by Stan.
 *
 * @tparam Evaluator an `Eigen::internal::evaluator<Xpr>` type
 */
template <typename Evaluator, typename = void>
struct is_stan_ad_product_hook : std::false_type {};
template <typename Evaluator>
struct is_stan_ad_product_hook<
    Evaluator, std::void_t<decltype(Evaluator::is_stan_ad_product_hook)>>
    : std::integral_constant<bool, Evaluator::is_stan_ad_product_hook> {};

/**
 * True if Eigen evaluates the expression `Xpr` with one of Stan's autodiff
 * product evaluators.
 *
 * @tparam Xpr an Eigen expression type
 */
template <typename Xpr>
inline constexpr bool is_stan_ad_product_hook_v
    = is_stan_ad_product_hook<Eigen::internal::evaluator<Xpr>>::value;

/**
 * Eigen's generic product evaluator for dense operands, reproduced from
 * Eigen/src/Core/ProductEvaluators.h (`product_evaluator<Product<Lhs, Rhs,
 * Options>, ProductTag, LhsShape, RhsShape>`). It is used for products of
 * scalars that Stan does not evaluate itself.
 *
 * @tparam Lhs left operand expression type
 * @tparam Rhs right operand expression type
 * @tparam Options Eigen product options
 * @tparam ProductTag Eigen product implementation type
 */
template <typename Lhs, typename Rhs, int Options, int ProductTag>
struct eigen_dense_product_evaluator
    : public Eigen::internal::evaluator<
          typename Eigen::Product<Lhs, Rhs, Options>::PlainObject> {
  using XprType = Eigen::Product<Lhs, Rhs, Options>;
  using PlainObject = typename XprType::PlainObject;
  using Base = Eigen::internal::evaluator<PlainObject>;
  enum { Flags = Base::Flags | Eigen::EvalBeforeNestingBit };

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE explicit eigen_dense_product_evaluator(
      const XprType& xpr)
      : m_result(xpr.rows(), xpr.cols()) {
    Eigen::internal::construct_at<Base>(this, m_result);
    Eigen::internal::generic_product_impl<Lhs, Rhs, Eigen::DenseShape,
                                          Eigen::DenseShape,
                                          ProductTag>::evalTo(m_result,
                                                              xpr.lhs(),
                                                              xpr.rhs());
  }

 protected:
  PlainObject m_result;
};

/**
 * `Assignment` of a product to a dense destination with the same autodiff
 * scalar. Eigen's own specializations for this case call
 * `generic_product_impl` directly and skip `product_evaluator`; for dense
 * operands this forwards to `eigen_matrix_product_impl` instead. Products
 * with non-dense operands (e.g. triangular views) run Eigen's body.
 *
 * @tparam Dst destination expression type
 * @tparam Lhs left operand expression type
 * @tparam Rhs right operand expression type
 * @tparam Options Eigen product options
 * @tparam Func assignment functor (`assign_op`, `add_assign_op`,
 * `sub_assign_op`)
 */
template <typename Dst, typename Lhs, typename Rhs, int Options, typename Func>
struct ad_product_assignment {
  using SrcXprType = Eigen::Product<Lhs, Rhs, Options>;
  using Scalar = typename Dst::Scalar;
  static constexpr bool dense_operands
      = std::is_same<typename Eigen::internal::evaluator_traits<Lhs>::Shape,
                     Eigen::DenseShape>::value
        && std::is_same<typename Eigen::internal::evaluator_traits<Rhs>::Shape,
                        Eigen::DenseShape>::value;
  static void run(Dst& dst, const SrcXprType& src, const Func& func) {
    if constexpr (dense_operands) {
      const auto result
          = eigen_matrix_product_impl<Lhs, Rhs>::run(src.lhs(), src.rhs());
      Eigen::internal::call_assignment_no_alias(dst, result, func);
    } else if constexpr (std::is_same<Func, Eigen::internal::assign_op<
                                                Scalar, Scalar>>::value) {
      if (dst.rows() != src.rows() || dst.cols() != src.cols()) {
        dst.resize(src.rows(), src.cols());
      }
      Eigen::internal::generic_product_impl<Lhs, Rhs>::evalTo(dst, src.lhs(),
                                                              src.rhs());
    } else if constexpr (std::is_same<Func, Eigen::internal::add_assign_op<
                                                Scalar, Scalar>>::value) {
      Eigen::internal::generic_product_impl<Lhs, Rhs>::addTo(dst, src.lhs(),
                                                             src.rhs());
    } else {
      Eigen::internal::generic_product_impl<Lhs, Rhs>::subTo(dst, src.lhs(),
                                                             src.rhs());
    }
  }
};

/**
 * Tag passed to Eigen's lazy (coefficient-based) product evaluator when it is
 * reached from one of Stan's per-tag specializations. Eigen's lazy evaluator
 * does not depend on its product tag, and no Eigen product has this tag, so
 * it resolves to Eigen's own evaluator.
 */
constexpr int eigen_lazy_product_passthrough_tag = -1;

template <typename Lhs, typename Rhs>
using eigen_lazy_product_evaluator = Eigen::internal::product_evaluator<
    Eigen::Product<Lhs, Rhs, Eigen::LazyProduct>,
    eigen_lazy_product_passthrough_tag, Eigen::DenseShape, Eigen::DenseShape>;

}  // namespace internal
}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

/**
 * Dense-by-dense products: Stan's evaluator for autodiff operands, Eigen's
 * generic evaluator otherwise. More specialized than Eigen's generic
 * `product_evaluator` (fixed shapes) and less specialized than Eigen's lazy
 * and `LazyCoeffBasedProductMode` evaluators, which are handled below.
 */
template <typename Lhs, typename Rhs, int Options, int ProductTag>
struct product_evaluator<Product<Lhs, Rhs, Options>, ProductTag, DenseShape,
                         DenseShape>
    : std::conditional_t<
          stan::math::internal::is_stan_ad_product_v<Lhs, Rhs>,
          stan::math::internal::ad_product_evaluator<Lhs, Rhs, Options>,
          stan::math::internal::eigen_dense_product_evaluator<Lhs, Rhs, Options,
                                                              ProductTag>> {
  using Base = std::conditional_t<
      stan::math::internal::is_stan_ad_product_v<Lhs, Rhs>,
      stan::math::internal::ad_product_evaluator<Lhs, Rhs, Options>,
      stan::math::internal::eigen_dense_product_evaluator<Lhs, Rhs, Options,
                                                          ProductTag>>;
  using Base::Base;
};

/**
 * Lazy (coefficient-based) products, one specialization per product tag so
 * each is more specialized than Eigen's `LazyProduct` evaluator. Eigen's
 * `LazyCoeffBasedProductMode` evaluator for default products derives from
 * the `CoeffBasedProductMode` one, so it is covered as well.
 */
#define STAN_EIGEN_LAZY_PRODUCT_EVALUATOR(TAG)                                 \
  template <typename Lhs, typename Rhs>                                        \
  struct product_evaluator<Product<Lhs, Rhs, LazyProduct>, TAG, DenseShape,    \
                           DenseShape>                                         \
      : std::conditional_t<                                                    \
            stan::math::internal::is_stan_ad_product_v<Lhs, Rhs>,              \
            stan::math::internal::ad_product_evaluator<Lhs, Rhs, LazyProduct>, \
            stan::math::internal::eigen_lazy_product_evaluator<Lhs, Rhs>> {    \
    using Base = std::conditional_t<                                           \
        stan::math::internal::is_stan_ad_product_v<Lhs, Rhs>,                  \
        stan::math::internal::ad_product_evaluator<Lhs, Rhs, LazyProduct>,     \
        stan::math::internal::eigen_lazy_product_evaluator<Lhs, Rhs>>;         \
    using Base::Base;                                                          \
  };

STAN_EIGEN_LAZY_PRODUCT_EVALUATOR(InnerProduct)
STAN_EIGEN_LAZY_PRODUCT_EVALUATOR(OuterProduct)
STAN_EIGEN_LAZY_PRODUCT_EVALUATOR(CoeffBasedProductMode)
STAN_EIGEN_LAZY_PRODUCT_EVALUATOR(LazyCoeffBasedProductMode)
STAN_EIGEN_LAZY_PRODUCT_EVALUATOR(GemvProduct)
STAN_EIGEN_LAZY_PRODUCT_EVALUATOR(GemmProduct)
#undef STAN_EIGEN_LAZY_PRODUCT_EVALUATOR

}  // namespace internal
}  // namespace Eigen

/**
 * Apply `X(TPL, LhsScalar, RhsScalar)` to every scalar pattern of a product
 * with an autodiff operand of type `AD` (and `std::complex<AD>`), where `S_`
 * stands for any other scalar. `TPL` and `TPL_S` are function-like macros
 * expanding to the extra template parameters (with a trailing comma) of the
 * patterns without and with `S_`. The `(AD, AD)`, `(complex<AD>,
 * complex<AD>)`, `(AD, complex<AD>)` and `(complex<AD>, AD)` patterns break
 * the ties between the patterns with `S_`.
 */
#define STAN_EIGEN_AD_SCALAR_PATTERNS(X, TPL, TPL_S, AD) \
  X(TPL_S, AD, S_)                                       \
  X(TPL_S, S_, AD)                                       \
  X(TPL, AD, AD)                                         \
  X(TPL_S, std::complex<AD>, S_)                         \
  X(TPL_S, S_, std::complex<AD>)                         \
  X(TPL, std::complex<AD>, std::complex<AD>)             \
  X(TPL, AD, std::complex<AD>)                           \
  X(TPL, std::complex<AD>, AD)

/**
 * Specialize Eigen's `binary_evaluator` for coefficient-wise products with
 * scalar pattern `(LS, RS)`, forwarding to `eigen_cwise_product_impl`. The
 * scalar template parameters of `binary_evaluator` stay defaulted, exactly as
 * in Eigen's generic specialization, so this one is more specialized only
 * through the functor type. Use inside `namespace Eigen::internal`.
 */
#define STAN_EIGEN_AD_CWISE_PRODUCT_EVALUATOR(TPL, LS, RS)                    \
  template <TPL() typename Lhs, typename Rhs>                                 \
  struct binary_evaluator<CwiseBinaryOp<scalar_product_op<LS, RS>, Lhs, Rhs>, \
                          IndexBased, IndexBased>                             \
      : stan::math::internal::ad_cwise_product_evaluator<                     \
            CwiseBinaryOp<scalar_product_op<LS, RS>, Lhs, Rhs>> {             \
    using Base = stan::math::internal::ad_cwise_product_evaluator<            \
        CwiseBinaryOp<scalar_product_op<LS, RS>, Lhs, Rhs>>;                  \
    using Base::Base;                                                         \
  };

/**
 * Specialize Eigen's `Assignment` of a product to a dense destination with
 * scalar `SCALAR` for `OP` (`assign_op`, `add_assign_op`, `sub_assign_op`).
 * The last template argument repeats Eigen's own condition token for token,
 * so this is more specialized only through the functor's scalar. Use inside
 * `namespace Eigen::internal`.
 */
#define STAN_EIGEN_AD_PRODUCT_ASSIGNMENT_OP(TPL, SCALAR, OP)            \
  template <TPL() typename DstXprType, typename Lhs, typename Rhs,      \
            int Options>                                                \
  struct Assignment<DstXprType, Product<Lhs, Rhs, Options>,             \
                    OP<SCALAR, SCALAR>, Dense2Dense,                    \
                    std::enable_if_t<(Options == DefaultProduct         \
                                      || Options == AliasFreeProduct)>> \
      : stan::math::internal::ad_product_assignment<                    \
            DstXprType, Lhs, Rhs, Options, OP<SCALAR, SCALAR>> {};

#define STAN_EIGEN_AD_PRODUCT_ASSIGNMENT(TPL, SCALAR)             \
  STAN_EIGEN_AD_PRODUCT_ASSIGNMENT_OP(TPL, SCALAR, assign_op)     \
  STAN_EIGEN_AD_PRODUCT_ASSIGNMENT_OP(TPL, SCALAR, add_assign_op) \
  STAN_EIGEN_AD_PRODUCT_ASSIGNMENT_OP(TPL, SCALAR, sub_assign_op)

#endif
