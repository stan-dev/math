#ifndef STAN_MATH_REV_FUN_INV_HPP
#define STAN_MATH_REV_FUN_INV_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/inv_square.hpp>

namespace stan {
namespace math {

namespace internal {
class inv_vari : public op_v_vari {
 public:
  explicit inv_vari(vari* avi) : op_v_vari(inv(avi->val_), avi) {}
  void chain() { avi_->adj_ -= adj_ * inv_square(avi_->val_); }
};

template <typename T>
class inv_matrix_vari : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int A_size_;
  double* Ad_;
  vari** variRefA_;
  vari** variRefInv_;

  /**
   * Constructor for exp_matrix_vari.
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
   */
  explicit inv_matrix_vari(const T& A)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        A_size_(A.size()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(A_size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)),
        variRefInv_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            A_size_)) {
    using Eigen::Map;
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_) = A.vi();
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Ad = A.val();
    Map<matrix_vi>(variRefInv_, A_rows_, A_cols_).array()
        = Ad.array().inverse()
                    .unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    Map<matrix_vi> RefInv(variRefInv_, A_rows_, A_cols_);
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj().array()
          -= RefInv.adj().array() * RefInv.val().array().square();
  }
};

}  // namespace internal


/**
 *
   \f[
   \mbox{inv}(x) =
   \begin{cases}
     \frac{1}{x} & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{inv}(x)}{\partial x} =
   \begin{cases}
     -\frac{1}{x^2} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 */
inline var inv(const var& a) { return var(new internal::inv_vari(a.vi_)); }

template <typename Container,
          require_container_st<is_container, is_var, Container>...>
inline auto inv(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_plain = plain_type_t<decltype(v)>;
        using T_ref = Eigen::Ref<const T_plain>;

        const T_ref& v_ref = v;
        auto* baseVari = new internal::inv_matrix_vari<T_ref>(v_ref);
        T_plain result(v_ref.rows(), v_ref.cols());
        result.vi() = Eigen::Map<matrix_vi>(baseVari->variRefInv_,
                                          v_ref.rows(), v_ref.cols());

        return result;
});
}

}  // namespace math
}  // namespace stan
#endif
