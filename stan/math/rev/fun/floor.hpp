#ifndef STAN_MATH_REV_FUN_FLOOR_HPP
#define STAN_MATH_REV_FUN_FLOOR_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
class floor_vari : public op_v_vari {
 public:
  explicit floor_vari(vari* avi) : op_v_vari(std::floor(avi->val_), avi) {}
  void chain() {
    if (unlikely(is_nan(avi_->val_))) {
      avi_->adj_ = NOT_A_NUMBER;
    }
  }
};

template <typename Container>
class floor_matrix_vari : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int A_size_;
  double* Ad_;
  vari** variRefA_;
  vari** variRefFloor_;

  /**
   * Constructor for floor_matrix_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * @tparam Container Type of Eigen expression/object
   * @param A Eigen expression/object
   */
  explicit floor_matrix_vari(const Container& A)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        A_size_(A.size()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(A_size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)),
        variRefFloor_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)) {
    using Eigen::Map;
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_) = A.vi();
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Ad = A.val();
    Map<matrix_vi>(variRefFloor_, A_rows_, A_cols_).array()
        = Ad.array().floor().unaryExpr(
            [](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj()
        = Ad.unaryExpr([](double x) { return std::isnan(x) ? x : 0.0; });
  }
};
}  // namespace internal

/**
 * Return the floor of the specified variable (cmath).
 *
 * The derivative of the floor function is defined and
 * zero everywhere but at integers, so we set these derivatives
 * to zero for convenience,
 *
 * \f$\frac{d}{dx} {\lfloor x \rfloor} = 0\f$.
 *
 * The floor function rounds down.  For double values, this is the largest
 * integral value that is not greater than the specified value.
 * Although this function is not differentiable because it is
 * discontinuous at integral values, its gradient is returned as
 * zero everywhere.
 *
   \f[
   \mbox{floor}(x) =
   \begin{cases}
     \lfloor x \rfloor & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{floor}(x)}{\partial x} =
   \begin{cases}
     0 & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Input variable.
 * @return Floor of the variable.
 */
inline var floor(const var& a) { return var(new internal::floor_vari(a.vi_)); }

/**
 * Return the floor of each variable in a container.
 *
 * @tparam Container Type of container
 * @param x Container of var
 * @return Floor of each variable in container.
 */
template <typename Container,
          require_container_st<is_container, is_var, Container>...>
inline auto floor(const Container& x) {
  return apply_vector_unary<Container>::apply(x, [](const auto& v) {
    using T_plain = plain_type_t<decltype(v)>;
    using T_ref = Eigen::Ref<const T_plain>;

    const T_ref& v_ref = v;
    auto* baseVari = new internal::floor_matrix_vari<T_ref>(v_ref);
    T_plain result(v_ref.rows(), v_ref.cols());
    result.vi() = Eigen::Map<matrix_vi>(baseVari->variRefFloor_, v_ref.rows(),
                                        v_ref.cols());

    return result;
  });
}
}  // namespace math
}  // namespace stan
#endif
