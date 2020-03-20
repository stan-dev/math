#ifndef STAN_MATH_REV_FUN_ACOS_HPP
#define STAN_MATH_REV_FUN_ACOS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
class acos_vari : public op_v_vari {
 public:
  explicit acos_vari(vari* avi) : op_v_vari(std::acos(avi->val_), avi) {}
  void chain() {
    avi_->adj_ -= adj_ / std::sqrt(1.0 - (avi_->val_ * avi_->val_));
  }
};

template <typename T>
class acos_matrix_vari : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int A_size_;
  double* Ad_;
  vari** variRefA_;
  vari** variRefAcos_;

  /**
   * Constructor for acos_matrix_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * @tparam T Type of Eigen expression/object
   * @param A Eigen expression/object
   */
  explicit acos_matrix_vari(const T& A)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        A_size_(A.size()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(A_size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)),
        variRefAcos_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            A_size_)) {
    using Eigen::Map;
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_) = A.vi();
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Ad = A.val();
    Map<matrix_vi>(variRefAcos_, A_rows_, A_cols_).array()
        = Ad.array().acos()
                    .unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    Map<matrix_vi> RefAcos(variRefAcos_, A_rows_, A_cols_);
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj().array()
          -= RefAcos.adj().array() * (1 - Ad.val().array().square()).rsqrt();
  }
};
}  // namespace internal

/**
 * Return the principal value of the arc cosine of a variable,
 * in radians (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \arccos x = \frac{-1}{\sqrt{1 - x^2}}\f$.
 *
 *
   \f[
   \mbox{acos}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \arccos(x) & \mbox{if } -1\leq x\leq 1 \\
     \textrm{NaN} & \mbox{if } x > 1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{acos}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \frac{\partial\, \arccos(x)}{\partial x} & \mbox{if } -1\leq x\leq 1 \\
     \textrm{NaN} & \mbox{if } x < -1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial \, \arccos(x)}{\partial x} = -\frac{1}{\sqrt{1-x^2}}
   \f]
 *
 * @param a Variable in range [-1, 1].
 * @return Arc cosine of variable, in radians.
 */
inline var acos(const var& a) { return var(new internal::acos_vari(a.vi_)); }

template <typename Container,
          require_container_st<is_container, is_var, Container>...>
inline auto acos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_plain = plain_type_t<decltype(v)>;
        using T_ref = Eigen::Ref<const T_plain>;

        const T_ref& v_ref = v;
        auto* baseVari = new internal::acos_matrix_vari<T_ref>(v_ref);
        T_plain result(v_ref.rows(), v_ref.cols());
        result.vi() = Eigen::Map<matrix_vi>(baseVari->variRefAcos_,
                                          v_ref.rows(), v_ref.cols());

        return result;
});
}

}  // namespace math
}  // namespace stan
#endif
