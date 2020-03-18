#ifndef STAN_MATH_REV_FUN_EXP_HPP
#define STAN_MATH_REV_FUN_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <cmath>
#include <iostream>

namespace stan {
namespace math {

namespace internal {
class exp_vari : public op_v_vari {
 public:
  explicit exp_vari(vari* avi) : op_v_vari(std::exp(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ * val_; }
};

template <typename T>
class exp_matrix_vari : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int A_size_;
  double* Ad_;
  vari** variRefA_;
  vari** variRefExp_;

  /**
   * Constructor for multiply_mat_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled by the second argument to
   * vari's constructor.
   *
   * @param A matrix
   * @param B matrix
   */
  explicit exp_matrix_vari(const T& A)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        A_size_(A.size()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(A_size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)),
        variRefExp_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            A_size_)) {
    using Eigen::Map;
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_) = A.vi();
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Ad = A.val();
    Map<matrix_vi>(variRefExp_, A_rows_, A_cols_).array()
        = Ad.array().exp().unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    Map<matrix_vi> RefExp(variRefExp_, A_rows_, A_cols_);
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj()
          += RefExp.val().cwiseProduct(RefExp.adj());
  }
};

}  // namespace internal



/*
template <int R, int C>
inline Eigen::Matrix<var, R, C> exp(const Eigen::Matrix<var, R, C>& A) {

  // Memory managed with the arena allocator.
  internal::exp_matrix_vari *baseVari = new internal::exp_matrix_vari(A);
  Eigen::Matrix<var, R, C> AB_v(A.rows(), A.cols());
  AB_v.vi() = Eigen::Map<matrix_vi>(baseVari->variRefAB_, A.rows(), A.cols());

  return AB_v;
}*/

template <typename Container,
          require_container_st<is_container, is_var, Container>...>
inline auto exp(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_vec = decltype(v);
        // Memory managed with the arena allocator.
        auto* baseVari = new internal::exp_matrix_vari<T_vec>(v);
        plain_type_t<T_vec> AB_v(v.rows(), v.cols());
        AB_v.vi() = Eigen::Map<matrix_vi>(baseVari->variRefExp_, v.rows(), v.cols());

        return AB_v;
});
}

/**
 * Return the exponentiation of the specified variable (cmath).
 *
   \f[
   \mbox{exp}(x) =
   \begin{cases}
     e^x & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{exp}(x)}{\partial x} =
   \begin{cases}
     e^x & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable to exponentiate.
 * @return Exponentiated variable.
 */
inline var exp(const var& a) { return var(new internal::exp_vari(a.vi_)); }

}  // namespace math
}  // namespace stan
#endif
