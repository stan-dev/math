#ifndef STAN_MATH_PRIM_META_EIGEN_PRODUCT_EVALUATOR_HOOKS_HPP
#define STAN_MATH_PRIM_META_EIGEN_PRODUCT_EVALUATOR_HOOKS_HPP

/**
 * @file
 * Specializations of Eigen's evaluators that compute products of matrices
 * with autodiff scalars with Stan's `multiply`, `elt_multiply`,
 * `scale_multiply` and `multiply_add`, so that `A * B` records one node on
 * the reverse-mode stack instead of one per scalar operation.
 *
 * Eigen's `product_evaluator` takes the operand scalar types as defaulted
 * template parameters, so a specialization that fixes them is ambiguous with
 * Eigen's generic one. The dense `product_evaluator` below therefore matches
 * every dense-by-dense product and handles non-autodiff scalars exactly as
 * Eigen's generic `product_evaluator` in Eigen/src/Core/ProductEvaluators.h
 * does. When upgrading Eigen, compare that evaluator with this one.
 *
 * Because the specializations match every dense product, this file is
 * included from stan/math/prim/fun/Eigen.hpp so that every translation unit
 * sees them. The functions they call are found by argument-dependent lookup
 * when a product is evaluated, so a translation unit that evaluates a `var`
 * product must include stan/math/rev.hpp.
 */

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_complex.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_var.hpp>
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
    = (is_var<typename Eigen::internal::traits<Lhs>::Scalar>::value
       || is_var<typename Eigen::internal::traits<Rhs>::Scalar>::value)
      && !is_complex<typename Eigen::internal::traits<Lhs>::Scalar>::value
      && !is_complex<typename Eigen::internal::traits<Rhs>::Scalar>::value;

/**
 * True for the expression Eigen builds for a scalar times a matrix, `s * X`:
 * a `CwiseBinaryOp` of `scalar_product_op` whose left operand is a
 * `CwiseNullaryOp<scalar_constant_op>` holding `s`
 * (`xpr.lhs().functor().m_other`) and whose right operand is `X` (`xpr.rhs()`).
 *
 * @tparam T an Eigen expression type
 */
template <typename T>
struct is_scalar_times_matrix : std::false_type {};
template <typename S1, typename S2, typename Plain, typename Mat>
struct is_scalar_times_matrix<Eigen::CwiseBinaryOp<
    Eigen::internal::scalar_product_op<S1, S2>,
    const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<S1>, Plain>,
    const Mat>> : std::true_type {};

/**
 * True for the expression Eigen builds for a matrix times a scalar, `X * s`:
 * a `CwiseBinaryOp` of `scalar_product_op` whose left operand is `X`
 * (`xpr.lhs()`) and whose right operand is a
 * `CwiseNullaryOp<scalar_constant_op>` holding `s`
 * (`xpr.rhs().functor().m_other`).
 *
 * @tparam T an Eigen expression type
 */
template <typename T>
struct is_matrix_times_scalar : std::false_type {};
template <typename S1, typename S2, typename Mat, typename Plain>
struct is_matrix_times_scalar<
    Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<S1, S2>, const Mat,
                         const Eigen::CwiseNullaryOp<
                             Eigen::internal::scalar_constant_op<S2>, Plain>>>
    : std::true_type {};

/**
 * True for an Eigen matrix product expression.
 *
 * @tparam T an Eigen expression type
 */
template <typename T>
struct is_eigen_product : std::false_type {};
template <typename Lhs, typename Rhs, int Options>
struct is_eigen_product<Eigen::Product<Lhs, Rhs, Options>> : std::true_type {};

}  // namespace internal
}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

/**
 * Dense-by-dense products. For operands with `var` scalars the product is
 * `scale_multiply(s, A, B)` when either operand is Eigen's `s * A` (Eigen
 * rewrites `s * (A * B)` as `(s * A) * B`) and `multiply(lhs, rhs)`
 * otherwise. For all other scalars this is Eigen's generic
 * `product_evaluator`. More specialized than Eigen's generic
 * `product_evaluator` (fixed shapes) and less specialized than its lazy and
 * `LazyCoeffBasedProductMode` evaluators, which are handled below.
 */
template <typename Lhs, typename Rhs, int Options, int ProductTag>
struct product_evaluator<Product<Lhs, Rhs, Options>, ProductTag, DenseShape,
                         DenseShape>
    : public evaluator<typename Product<Lhs, Rhs, Options>::PlainObject> {
  using XprType = Product<Lhs, Rhs, Options>;
  using PlainObject = typename XprType::PlainObject;
  using Base = evaluator<PlainObject>;
  enum { Flags = Base::Flags | EvalBeforeNestingBit };

  explicit product_evaluator(const XprType& xpr)
      : m_result(xpr.rows(), xpr.cols()) {
    using stan::math::internal::is_matrix_times_scalar;
    using stan::math::internal::is_scalar_times_matrix;
    if constexpr (!stan::math::internal::is_stan_ad_product_v<Lhs, Rhs>) {
      generic_product_impl<Lhs, Rhs, DenseShape, DenseShape,
                           ProductTag>::evalTo(m_result, xpr.lhs(), xpr.rhs());
    } else if constexpr (is_scalar_times_matrix<Lhs>::value) {
      // (s * A) * B
      m_result = scale_multiply(xpr.lhs().lhs().functor().m_other,
                                xpr.lhs().rhs(), xpr.rhs());
    } else if constexpr (is_matrix_times_scalar<Lhs>::value) {
      // (A * s) * B
      m_result = scale_multiply(xpr.lhs().rhs().functor().m_other,
                                xpr.lhs().lhs(), xpr.rhs());
    } else if constexpr (is_scalar_times_matrix<Rhs>::value) {
      // A * (s * B)
      m_result = scale_multiply(xpr.rhs().lhs().functor().m_other, xpr.lhs(),
                                xpr.rhs().rhs());
    } else if constexpr (is_matrix_times_scalar<Rhs>::value) {
      // A * (B * s)
      m_result = scale_multiply(xpr.rhs().rhs().functor().m_other, xpr.lhs(),
                                xpr.rhs().lhs());
    } else if constexpr (stan::is_stan_scalar<decltype(
                             multiply(xpr.lhs(), xpr.rhs()))>::value) {
      // a row vector times a column vector
      m_result.coeffRef(0, 0) = multiply(xpr.lhs(), xpr.rhs());
    } else {
      m_result = multiply(xpr.lhs(), xpr.rhs());
    }
    // Base points into m_result's storage, so it is constructed only once
    // m_result is filled: assigning a returned matrix can move a new buffer
    // into m_result.
    construct_at<Base>(this, m_result);
  }

 protected:
  PlainObject m_result;
};

/**
 * Lazy (coefficient-based) products, one specialization per product tag so
 * each is more specialized than Eigen's `LazyProduct` evaluator. For `var`
 * operands they use the dense evaluator above (the product tag is unused
 * there); for other scalars Eigen's lazy evaluator, reached with tag -1,
 * which no product has and which Eigen's lazy evaluator ignores. Eigen's
 * `LazyCoeffBasedProductMode` evaluator for default products derives from
 * the `CoeffBasedProductMode` specialization, so it is covered as well.
 */
#define STAN_EIGEN_LAZY_PRODUCT_EVALUATOR(TAG)                                \
  template <typename Lhs, typename Rhs>                                       \
  struct product_evaluator<Product<Lhs, Rhs, LazyProduct>, TAG, DenseShape,   \
                           DenseShape>                                        \
      : std::conditional_t<                                                   \
            stan::math::internal::is_stan_ad_product_v<Lhs, Rhs>,             \
            product_evaluator<Product<Lhs, Rhs, DefaultProduct>, GemmProduct, \
                              DenseShape, DenseShape>,                        \
            product_evaluator<Product<Lhs, Rhs, LazyProduct>, -1, DenseShape, \
                              DenseShape>> {                                  \
    using XprType = Product<Lhs, Rhs, LazyProduct>;                           \
    using Base = std::conditional_t<                                          \
        stan::math::internal::is_stan_ad_product_v<Lhs, Rhs>,                 \
        product_evaluator<Product<Lhs, Rhs, DefaultProduct>, GemmProduct,     \
                          DenseShape, DenseShape>,                            \
        product_evaluator<Product<Lhs, Rhs, LazyProduct>, -1, DenseShape,     \
                          DenseShape>>;                                       \
    explicit product_evaluator(const XprType& xpr) : Base(base_xpr(xpr)) {}   \
                                                                              \
   private:                                                                   \
    static auto base_xpr(const XprType& xpr) {                                \
      if constexpr (stan::math::internal::is_stan_ad_product_v<Lhs, Rhs>) {   \
        return Product<Lhs, Rhs, DefaultProduct>(xpr.lhs(), xpr.rhs());       \
      } else {                                                                \
        return xpr;                                                           \
      }                                                                       \
    }                                                                         \
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
 * Specialize Eigen's `binary_evaluator` for the coefficient-wise product
 * `lhs .* rhs` with scalar pattern `(LS, RS)`. The product is
 * `scale_multiply(s, A, B)` for `s * (A * B)` or `(A * B) * s`,
 * `multiply(s, X)` for `s * X` or `X * s`, and `elt_multiply(lhs, rhs)`
 * otherwise. `TPL` is a
 * function-like macro expanding to extra template parameters with a
 * trailing comma. The scalar template parameters of `binary_evaluator` stay
 * defaulted, as in Eigen's generic specialization, so this one is more
 * specialized only through the functor. Use inside `namespace
 * Eigen::internal`.
 */
#define STAN_EIGEN_AD_CWISE_PRODUCT_EVALUATOR(TPL, LS, RS)                    \
  template <TPL() typename Lhs, typename Rhs>                                 \
  struct binary_evaluator<CwiseBinaryOp<scalar_product_op<LS, RS>, Lhs, Rhs>, \
                          IndexBased, IndexBased>                             \
      : public evaluator<typename CwiseBinaryOp<scalar_product_op<LS, RS>,    \
                                                Lhs, Rhs>::PlainObject> {     \
    using XprType = CwiseBinaryOp<scalar_product_op<LS, RS>, Lhs, Rhs>;       \
    using PlainObject = typename XprType::PlainObject;                        \
    using Base = evaluator<PlainObject>;                                      \
    enum { Flags = Base::Flags | EvalBeforeNestingBit };                      \
                                                                              \
    explicit binary_evaluator(const XprType& xpr)                             \
        : m_result(xpr.rows(), xpr.cols()) {                                  \
      using stan::math::internal::is_eigen_product;                           \
      using stan::math::internal::is_matrix_times_scalar;                     \
      using stan::math::internal::is_scalar_times_matrix;                     \
      /* Stan's functions take matrices, so array operands are viewed as */   \
      /* matrices (`.matrix().derived()` is the operand itself for a */       \
      /* matrix) and the result is viewed back as an array */                 \
      const auto product = [&xpr]() {                                         \
        if constexpr (is_scalar_times_matrix<XprType>::value) {               \
          if constexpr (is_eigen_product<std::decay_t<Rhs>>::value) {         \
            /* s * (A * B) */                                                 \
            return scale_multiply(xpr.lhs().functor().m_other,                \
                                  xpr.rhs().lhs(), xpr.rhs().rhs());          \
          } else {                                                            \
            /* s * X */                                                       \
            return multiply(xpr.lhs().functor().m_other,                      \
                            xpr.rhs().matrix().derived());                    \
          }                                                                   \
        } else if constexpr (is_matrix_times_scalar<XprType>::value) {        \
          if constexpr (is_eigen_product<std::decay_t<Lhs>>::value) {         \
            /* (A * B) * s */                                                 \
            return scale_multiply(xpr.rhs().functor().m_other,                \
                                  xpr.lhs().lhs(), xpr.lhs().rhs());          \
          } else {                                                            \
            /* X * s */                                                       \
            return multiply(xpr.rhs().functor().m_other,                      \
                            xpr.lhs().matrix().derived());                    \
          }                                                                   \
        } else {                                                              \
          return elt_multiply(xpr.lhs().matrix().derived(),                   \
                              xpr.rhs().matrix().derived());                  \
        }                                                                     \
      }();                                                                    \
      if constexpr (std::is_base_of<ArrayBase<PlainObject>,                   \
                                    PlainObject>::value) {                    \
        m_result = product.array();                                           \
      } else {                                                                \
        m_result = product;                                                   \
      }                                                                       \
      /* Base points into m_result's storage, so it is constructed only */    \
      /* once m_result is filled */                                           \
      construct_at<Base>(this, m_result);                                     \
    }                                                                         \
                                                                              \
   protected:                                                                 \
    PlainObject m_result;                                                     \
  };

/**
 * Specialize Eigen's `Assignment` of a product to a destination with scalar
 * `SCALAR` for `OP` (`assign_op`, `add_assign_op` or `sub_assign_op`).
 * Eigen's own specializations for this case compute the product directly
 * instead of through `product_evaluator`. For dense operands, `dst = A * B`
 * uses `multiply` (or `scale_multiply` for a scaled operand), `dst += A * B`
 * uses `multiply_add(A, B, dst)` and `dst -= A * B` uses `subtract(dst,
 * multiply(A, B))`; other shapes (e.g. triangular views) run Eigen's body.
 * The last template argument repeats Eigen's own condition token for token,
 * so this is more specialized only through the functor's scalar. Use inside
 * `namespace Eigen::internal`.
 */
#define STAN_EIGEN_AD_PRODUCT_ASSIGNMENT_OP(TPL, SCALAR, OP)                  \
  template <TPL() typename DstXprType, typename Lhs, typename Rhs,            \
            int Options>                                                      \
  struct Assignment<DstXprType, Product<Lhs, Rhs, Options>,                   \
                    OP<SCALAR, SCALAR>, Dense2Dense,                          \
                    std::enable_if_t<(Options == DefaultProduct               \
                                      || Options == AliasFreeProduct)>> {     \
    using SrcXprType = Product<Lhs, Rhs, Options>;                            \
    static void run(DstXprType& dst, const SrcXprType& src,                   \
                    const OP<SCALAR, SCALAR>& func) {                         \
      using stan::math::internal::is_matrix_times_scalar;                     \
      using stan::math::internal::is_scalar_times_matrix;                     \
      constexpr bool is_add                                                   \
          = std::is_same<OP<SCALAR, SCALAR>,                                  \
                         add_assign_op<SCALAR, SCALAR>>::value;               \
      constexpr bool is_sub                                                   \
          = std::is_same<OP<SCALAR, SCALAR>,                                  \
                         sub_assign_op<SCALAR, SCALAR>>::value;               \
      constexpr bool is_scaled = is_scalar_times_matrix<Lhs>::value           \
                                 || is_matrix_times_scalar<Lhs>::value        \
                                 || is_scalar_times_matrix<Rhs>::value        \
                                 || is_matrix_times_scalar<Rhs>::value;       \
      if constexpr (!std::is_same<typename evaluator_traits<Lhs>::Shape,      \
                                  DenseShape>::value                          \
                    || !std::is_same<typename evaluator_traits<Rhs>::Shape,   \
                                     DenseShape>::value) {                    \
        if constexpr (is_add) {                                               \
          generic_product_impl<Lhs, Rhs>::addTo(dst, src.lhs(), src.rhs());   \
        } else if constexpr (is_sub) {                                        \
          generic_product_impl<Lhs, Rhs>::subTo(dst, src.lhs(), src.rhs());   \
        } else {                                                              \
          if (dst.rows() != src.rows() || dst.cols() != src.cols()) {         \
            dst.resize(src.rows(), src.cols());                               \
          }                                                                   \
          generic_product_impl<Lhs, Rhs>::evalTo(dst, src.lhs(), src.rhs());  \
        }                                                                     \
      } else if constexpr (stan::is_stan_scalar<decltype(                     \
                               multiply(src.lhs(), src.rhs()))>::value) {     \
        /* a row vector times a column vector */                              \
        if (dst.rows() != 1 || dst.cols() != 1) {                             \
          dst.resize(1, 1);                                                   \
        }                                                                     \
        func.assignCoeff(dst.coeffRef(0, 0), multiply(src.lhs(), src.rhs())); \
      } else if constexpr (is_add && !is_scaled) {                            \
        dst = multiply_add(src.lhs(), src.rhs(), dst);                        \
      } else if constexpr (is_sub && !is_scaled) {                            \
        dst = subtract(dst, multiply(src.lhs(), src.rhs()));                  \
      } else {                                                                \
        const auto product = [&src]() {                                       \
          if constexpr (is_scalar_times_matrix<Lhs>::value) {                 \
            /* (s * A) * B */                                                 \
            return scale_multiply(src.lhs().lhs().functor().m_other,          \
                                  src.lhs().rhs(), src.rhs());                \
          } else if constexpr (is_matrix_times_scalar<Lhs>::value) {          \
            /* (A * s) * B */                                                 \
            return scale_multiply(src.lhs().rhs().functor().m_other,          \
                                  src.lhs().lhs(), src.rhs());                \
          } else if constexpr (is_scalar_times_matrix<Rhs>::value) {          \
            /* A * (s * B) */                                                 \
            return scale_multiply(src.rhs().lhs().functor().m_other,          \
                                  src.lhs(), src.rhs().rhs());                \
          } else if constexpr (is_matrix_times_scalar<Rhs>::value) {          \
            /* A * (B * s) */                                                 \
            return scale_multiply(src.rhs().rhs().functor().m_other,          \
                                  src.lhs(), src.rhs().lhs());                \
          } else {                                                            \
            return multiply(src.lhs(), src.rhs());                            \
          }                                                                   \
        }();                                                                  \
        if constexpr (is_add) {                                               \
          dst += product;                                                     \
        } else if constexpr (is_sub) {                                        \
          dst -= product;                                                     \
        } else {                                                              \
          dst = product;                                                      \
        }                                                                     \
      }                                                                       \
    }                                                                         \
  };

#define STAN_EIGEN_AD_PRODUCT_ASSIGNMENT(TPL, SCALAR)             \
  STAN_EIGEN_AD_PRODUCT_ASSIGNMENT_OP(TPL, SCALAR, assign_op)     \
  STAN_EIGEN_AD_PRODUCT_ASSIGNMENT_OP(TPL, SCALAR, add_assign_op) \
  STAN_EIGEN_AD_PRODUCT_ASSIGNMENT_OP(TPL, SCALAR, sub_assign_op)

#endif
