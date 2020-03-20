#ifndef STAN_MATH_REV_FUN_ASIN_HPP
#define STAN_MATH_REV_FUN_ASIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
class asin_vari : public op_v_vari {
 public:
  explicit asin_vari(vari* avi) : op_v_vari(std::asin(avi->val_), avi) {}
  void chain() {
    avi_->adj_ += adj_ / std::sqrt(1.0 - (avi_->val_ * avi_->val_));
  }
};

template <typename Container>
class asin_matrix_vari : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int A_size_;
  double* Ad_;
  vari** variRefA_;
  vari** variRefAsin_;

  /**
   * Constructor for asin_matrix_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * @tparam Container Type of Eigen expression/object
   * @param A Eigen expression/object
   */
  explicit asin_matrix_vari(const Container& A)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        A_size_(A.size()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(A_size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)),
        variRefAsin_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)) {
    using Eigen::Map;
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_) = A.vi();
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Ad = A.val();
    Map<matrix_vi>(variRefAsin_, A_rows_, A_cols_).array()
        = Ad.array().asin().unaryExpr(
            [](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    Map<matrix_vi> RefAsin(variRefAsin_, A_rows_, A_cols_);
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj().array()
        += RefAsin.adj().array() * (1 - Ad.array().square()).rsqrt();
  }
};
}  // namespace internal

/**
 * Return the principal value of the arc sine, in radians, of the
 * specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \arcsin x = \frac{1}{\sqrt{1 - x^2}}\f$.
 *
 *
   \f[
   \mbox{asin}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \arcsin(x) & \mbox{if } -1\leq x\leq 1 \\
     \textrm{NaN} & \mbox{if } x > 1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{asin}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \frac{\partial\, \arcsin(x)}{\partial x} & \mbox{if } -1\leq x\leq 1 \\
     \textrm{NaN} & \mbox{if } x < -1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial \, \arcsin(x)}{\partial x} = \frac{1}{\sqrt{1-x^2}}
   \f]
 *
 * @param a Variable in range [-1, 1].
 * @return Arc sine of variable, in radians.
 */
inline var asin(const var& a) { return var(new internal::asin_vari(a.vi_)); }

/**
 * Return the arcsine of each variable in a container.
 *
 * @tparam Container Type of container
 * @param x Container of var
 * @return Arcsine of each variable in container.
 */
template <typename Container,
          require_container_st<is_container, is_var, Container>...>
inline auto asin(const Container& x) {
  return apply_vector_unary<Container>::apply(x, [](const auto& v) {
    using T_plain = plain_type_t<decltype(v)>;
    using T_ref = Eigen::Ref<const T_plain>;

    const T_ref& v_ref = v;
    auto* baseVari = new internal::asin_matrix_vari<T_ref>(v_ref);
    T_plain result(v_ref.rows(), v_ref.cols());
    result.vi() = Eigen::Map<matrix_vi>(baseVari->variRefAsin_, v_ref.rows(),
                                        v_ref.cols());

    return result;
  });
}
}  // namespace math
}  // namespace stan
#endif
