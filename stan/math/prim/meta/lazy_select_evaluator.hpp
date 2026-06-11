#ifndef STAN_MATH_PRIM_META_LAZY_SELECT_EVALUATOR_HPP
#define STAN_MATH_PRIM_META_LAZY_SELECT_EVALUATOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Evaluator for `.select()` expressions that evaluates only the branch
 * chosen by the condition, restoring the Eigen 3.x `Select` semantics.
 *
 * Eigen 5 implements `select` as a `CwiseTernaryOp` whose evaluator
 * computes the then-branch, else-branch, and condition coefficients
 * eagerly before applying the `cond == 0 ? b : a` functor. For plain
 * arithmetic scalars the discarded branch is dead arithmetic, but for
 * reverse-mode autodiff scalars evaluating a branch has a side effect:
 * it allocates varis on the autodiff stack. `chain()` runs on every vari
 * during the reverse pass, including orphans whose adjoint is zero, so an
 * orphan whose chain rule multiplies its adjoint by a value-derived
 * quantity (e.g. `exp`: `adj * exp(val)` with `val = +INF`) injects
 * `0 * INF = NaN` into the adjoint of a live input. Distribution code
 * throughout `stan/math/prim/prob` relies on `.select()` to guard exactly
 * such numerically dangerous expressions.
 *
 * This base class provides the lazy coefficient access; the
 * `Eigen::internal::ternary_evaluator` partial specializations for
 * `scalar_boolean_select_op` on `var` (rev/core/Eigen_NumTraits.hpp) and
 * `fvar<T>` (fwd/fun/Eigen_NumTraits.hpp) derive from it. Specializing
 * `ternary_evaluator` also intercepts Eigen's fused
 * `(a < b).select(c, d)` evaluator, which derives from
 * `ternary_evaluator` after rebuilding the expression.
 *
 * The packet (vectorized) path needs no override: for autodiff scalars
 * `functor_traits<scalar_boolean_select_op>::PacketAccess` is false, and
 * `Flags` below additionally masks out `PacketAccessBit`.
 *
 * @tparam TernaryOp the select functor,
 * `Eigen::internal::scalar_boolean_select_op<Scalar, Scalar, CondScalar>`
 * @tparam Arg1 type of the then-branch expression
 * @tparam Arg2 type of the else-branch expression
 * @tparam Arg3 type of the condition expression
 */
template <typename TernaryOp, typename Arg1, typename Arg2, typename Arg3>
struct lazy_select_evaluator
    : Eigen::internal::evaluator_base<
          Eigen::CwiseTernaryOp<TernaryOp, Arg1, Arg2, Arg3>> {
  using XprType = Eigen::CwiseTernaryOp<TernaryOp, Arg1, Arg2, Arg3>;
  using CoeffReturnType = typename XprType::CoeffReturnType;
  using CondScalar = typename Arg3::Scalar;
  enum {
    CoeffReadCost
    = static_cast<int>(Eigen::internal::evaluator<Arg1>::CoeffReadCost)
      + static_cast<int>(Eigen::internal::evaluator<Arg2>::CoeffReadCost)
      + static_cast<int>(Eigen::internal::evaluator<Arg3>::CoeffReadCost)
      + static_cast<int>(Eigen::internal::functor_traits<TernaryOp>::Cost),

    Arg1Flags = Eigen::internal::evaluator<Arg1>::Flags,
    Arg2Flags = Eigen::internal::evaluator<Arg2>::Flags,
    Arg3Flags = Eigen::internal::evaluator<Arg3>::Flags,
    StorageOrdersAgree
    = (static_cast<int>(Arg1Flags) & Eigen::RowMajorBit)
          == (static_cast<int>(Arg2Flags) & Eigen::RowMajorBit)
      && (static_cast<int>(Arg1Flags) & Eigen::RowMajorBit)
             == (static_cast<int>(Arg3Flags) & Eigen::RowMajorBit),
    Flags0 = (static_cast<int>(Arg1Flags) | static_cast<int>(Arg2Flags)
              | static_cast<int>(Arg3Flags))
             & (Eigen::HereditaryBits
                | (static_cast<int>(Arg1Flags) & static_cast<int>(Arg2Flags)
                   & static_cast<int>(Arg3Flags)
                   & (StorageOrdersAgree ? Eigen::LinearAccessBit : 0))),
    Flags = (Flags0 & ~Eigen::RowMajorBit) | (Arg1Flags & Eigen::RowMajorBit),
    Alignment = Eigen::internal::plain_enum_min(
        Eigen::internal::plain_enum_min(
            Eigen::internal::evaluator<Arg1>::Alignment,
            Eigen::internal::evaluator<Arg2>::Alignment),
        Eigen::internal::evaluator<Arg3>::Alignment)
  };

  EIGEN_DEVICE_FUNC explicit lazy_select_evaluator(const XprType& xpr)
      : m_d(xpr) {}

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE CoeffReturnType
  coeff(Eigen::Index row, Eigen::Index col) const {
    return m_d.arg3Impl.coeff(row, col) == CondScalar(0)
               ? m_d.arg2Impl.coeff(row, col)
               : m_d.arg1Impl.coeff(row, col);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE CoeffReturnType
  coeff(Eigen::Index index) const {
    return m_d.arg3Impl.coeff(index) == CondScalar(0)
               ? m_d.arg2Impl.coeff(index)
               : m_d.arg1Impl.coeff(index);
  }

 protected:
  struct Data {
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE explicit Data(const XprType& xpr)
        : op(xpr.functor()),
          arg1Impl(xpr.arg1()),
          arg2Impl(xpr.arg2()),
          arg3Impl(xpr.arg3()) {}
    TernaryOp op;
    Eigen::internal::evaluator<Arg1> arg1Impl;
    Eigen::internal::evaluator<Arg2> arg2Impl;
    Eigen::internal::evaluator<Arg3> arg3Impl;
  };

  Data m_d;
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
