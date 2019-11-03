#ifndef STAN_MATH_REV_MAT_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/LDLT_alloc.hpp>
#include <stan/math/rev/mat/fun/LDLT_factor.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {
template <typename T2, int R2, int C2, typename T3, int R3, int C3>
class trace_inv_quad_form_ldlt_impl : public chainable_alloc {
 protected:
  inline void initializeB(const Eigen::Matrix<var, R3, C3> &B, bool haveD) {
    matrix_d Bd = B.val();
    variB_ = B.vi();
    AinvB_ = ldlt_.solve(Bd);
    if (haveD) {
      C_.noalias() = Bd.transpose() * AinvB_;
    } else {
      value_ = (Bd.transpose() * AinvB_).trace();
    }
  }
  inline void initializeB(const Eigen::Matrix<double, R3, C3> &B, bool haveD) {
    AinvB_ = ldlt_.solve(B);
    if (haveD) {
      C_.noalias() = B.transpose() * AinvB_;
    } else {
      value_ = (B.transpose() * AinvB_).trace();
    }
  }

  template <int R1, int C1>
  inline void initializeD(const Eigen::Matrix<var, R1, C1> &D) {
    D_ = D.val();
    variD_ = D.vi();
  }
  template <int R1, int C1>
  inline void initializeD(const Eigen::Matrix<double, R1, C1> &D) {
    D_ = D;
  }

 public:
  template <typename T1, int R1, int C1>
  trace_inv_quad_form_ldlt_impl(const Eigen::Matrix<T1, R1, C1> &D,
                                const LDLT_factor<T2, R2, C2> &A,
                                const Eigen::Matrix<T3, R3, C3> &B)
      : Dtype_(stan::is_var<T1>::value), ldlt_(A) {
    initializeB(B, true);
    initializeD(D);

    value_ = (D_ * C_).trace();
  }

  trace_inv_quad_form_ldlt_impl(const LDLT_factor<T2, R2, C2> &A,
                                const Eigen::Matrix<T3, R3, C3> &B)
      : Dtype_(2), ldlt_(A) {
    initializeB(B, false);
  }

  const int Dtype_;  // 0 = double, 1 = var, 2 = missing
  LDLT_factor<T2, R2, C2> ldlt_;
  matrix_d D_;
  matrix_vi variD_;
  matrix_vi variB_;
  matrix_d AinvB_;
  matrix_d C_;
  double value_;
};

template <typename T2, int R2, int C2, typename T3, int R3, int C3>
class trace_inv_quad_form_ldlt_vari : public vari {
 protected:
  static inline void chainA(
      double adj,
      trace_inv_quad_form_ldlt_impl<double, R2, C2, T3, R3, C3> *impl) {}
  static inline void chainB(
      double adj,
      trace_inv_quad_form_ldlt_impl<T2, R2, C2, double, R3, C3> *impl) {}

  static inline void chainA(
      double adj,
      trace_inv_quad_form_ldlt_impl<var, R2, C2, T3, R3, C3> *impl) {
    Eigen::Matrix<double, R2, C2> aA;

    if (impl->Dtype_ != 2) {
      aA.noalias()
          = -adj
            * (impl->AinvB_ * impl->D_.transpose() * impl->AinvB_.transpose());
    } else {
      aA.noalias() = -adj * (impl->AinvB_ * impl->AinvB_.transpose());
    }

    impl->ldlt_.alloc_->variA_.adj() += aA;
  }
  static inline void chainB(
      double adj,
      trace_inv_quad_form_ldlt_impl<T2, R2, C2, var, R3, C3> *impl) {
    matrix_d aB;

    if (impl->Dtype_ != 2) {
      aB.noalias() = adj * impl->AinvB_ * (impl->D_ + impl->D_.transpose());
    } else {
      aB.noalias() = 2 * adj * impl->AinvB_;
    }

    impl->variB_.adj() += aB;
  }

 public:
  explicit trace_inv_quad_form_ldlt_vari(
      trace_inv_quad_form_ldlt_impl<T2, R2, C2, T3, R3, C3> *impl)
      : vari(impl->value_), impl_(impl) {}

  virtual void chain() {
    // F = trace(D * B' * inv(A) * B)
    // aA = -aF * inv(A') * B * D' * B' * inv(A')
    // aB = aF*(inv(A) * B * D + inv(A') * B * D')
    // aD = aF*(B' * inv(A) * B)
    chainA(adj_, impl_);

    chainB(adj_, impl_);

    if (impl_->Dtype_ == 1) {
      impl_->variD_.adj() += adj_ * impl_->C_;
    }
  }

  trace_inv_quad_form_ldlt_impl<T2, R2, C2, T3, R3, C3> *impl_;
};

}  // namespace internal

/**
 * Compute the trace of an inverse quadratic form.  I.E., this computes
 *       trace(B^T A^-1 B)
 * where the LDLT_factor of A is provided.
 **/
template <typename T2, int R2, int C2, typename T3, int R3, int C3,
          typename = require_any_var_t<T2, T3>>
inline return_type_t<T2, T3> trace_inv_quad_form_ldlt(
    const LDLT_factor<T2, R2, C2> &A, const Eigen::Matrix<T3, R3, C3> &B) {
  check_multiplicable("trace_inv_quad_form_ldlt", "A", A, "B", B);

  internal::trace_inv_quad_form_ldlt_impl<T2, R2, C2, T3, R3, C3> *impl_
      = new internal::trace_inv_quad_form_ldlt_impl<T2, R2, C2, T3, R3, C3>(A,
                                                                            B);

  return var(
      new internal::trace_inv_quad_form_ldlt_vari<T2, R2, C2, T3, R3, C3>(
          impl_));
}

}  // namespace math
}  // namespace stan
#endif
