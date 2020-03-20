#ifndef STAN_MATH_REV_FUN_FABS_HPP
#define STAN_MATH_REV_FUN_FABS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

namespace internal {

template <typename Container>
class fabs_matrix_vari : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int A_size_;
  double* Ad_;
  vari** variRefA_;
  vari** variRefFabs_;

  /**
   * Constructor for fabs_matrix_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * @tparam Container Type of Eigen expression/object
   * @param A Eigen expression/object
   */
  explicit fabs_matrix_vari(const Container& A)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        A_size_(A.size()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(A_size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)),
        variRefFabs_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)) {
    using Eigen::Map;
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_) = A.vi();
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Ad = A.val();
    Map<matrix_vi>(variRefFabs_, A_rows_, A_cols_).array()
        = Ad.array().abs().unaryExpr(
            [](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Map<matrix_vi> RefFabs(variRefFabs_, A_rows_, A_cols_);
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj()
        = Ad.array().isNaN().select(
            NOT_A_NUMBER,
            (Ad.array() < 0)
                .select(-RefFabs.adj().array(), RefFabs.adj().array()));
  }
};
}  // namespace internal
/**
 * Return the absolute value of the variable (cmath).
 *
 * Choosing an arbitrary value at the non-differentiable point 0,
 *
 * \f$\frac{d}{dx}|x| = \mbox{sgn}(x)\f$.
 *
 * where \f$\mbox{sgn}(x)\f$ is the signum function, taking values
 * -1 if \f$x < 0\f$, 0 if \f$x == 0\f$, and 1 if \f$x == 1\f$.
 *
 * The function <code>abs()</code> provides the same behavior, with
 * <code>abs()</code> defined in stdlib.h and <code>fabs()</code>
 * defined in <code>cmath</code>.
 * The derivative is 0 if the input is 0.
 *
 * Returns std::numeric_limits<double>::quiet_NaN() for NaN inputs.
 *
 *
   \f[
   \mbox{fabs}(x) =
   \begin{cases}
     |x| & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fabs}(x)}{\partial x} =
   \begin{cases}
     -1 & \mbox{if } x < 0 \\
     0 & \mbox{if } x = 0 \\
     1 & \mbox{if } x > 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Input variable.
 * @return Absolute value of variable.
 */
inline var fabs(const var& a) {
  if (a.val() > 0.0) {
    return a;
  } else if (a.val() < 0.0) {
    return var(new internal::neg_vari(a.vi_));
  } else if (a.val() == 0) {
    return var(new vari(0));
  } else {
    return var(new precomp_v_vari(NOT_A_NUMBER, a.vi_, NOT_A_NUMBER));
  }
}

/**
 * Return the absolute value of each variable in a container.
 *
 * @tparam Container Type of container
 * @param x Container of var
 * @return Absolute value of each variable in container.
 */
template <typename Container,
          require_container_st<is_container, is_var, Container>...>
inline auto fabs(const Container& x) {
  return apply_vector_unary<Container>::apply(x, [](const auto& v) {
    using T_plain = plain_type_t<decltype(v)>;
    using T_ref = Eigen::Ref<const T_plain>;

    const T_ref& v_ref = v;
    auto* baseVari = new internal::fabs_matrix_vari<T_ref>(v_ref);
    T_plain result(v_ref.rows(), v_ref.cols());
    result.vi() = Eigen::Map<matrix_vi>(baseVari->variRefFabs_, v_ref.rows(),
                                        v_ref.cols());

    return result;
  });
}

}  // namespace math
}  // namespace stan
#endif
