#ifndef STAN_MATH_REV_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/adj_jac_apply.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename T>
class positive_ordered_constrain_op {
 public:
  adj_arg_t<T> exp_x_;
  explicit positive_ordered_constrain_op(const T& x) : exp_x_(setup_adj_arg<T>(x.size())) {}
  /**
   * Return an increasing positive ordered vector derived from the specified
   * free vector.  The returned constrained vector will have the
   * same dimensionality as the specified free vector.
   * @tparam ColVec An eigen type with compile time dynamic rows and 1 column
   * @param x Free vector of scalars
   * @return Positive, increasing ordered vector
   */
  template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
  inline Eigen::VectorXd forward_pass(ColVec&& x) {
    ref_type_t<ColVec> x_ref(std::forward<ColVec>(x));
    using std::exp;
    Eigen::Matrix<double, Eigen::Dynamic, 1> y(x_ref.size());
    if (x_ref.size() == 0) {
      return y;
    }
    exp_x_ = exp(x_ref.array());
    y[0] = exp_x_(0);
    for (int n = 1; n < x.size(); ++n) {
      y[n] = y[n - 1] + exp_x_(n);
    }
    return y;
  }

  /**
   * Compute the result of multiply the transpose of the adjoint vector times
   * the Jacobian of the positive_ordered_constrain operator.
   *
   * @tparam ColVec An eigen type with compile time dynamic rows and 1 column
   * @param adj Eigen::VectorXd of adjoints at the output of the softmax
   * @return Eigen::VectorXd of adjoints propagated through softmax operation
   */
  template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
  inline auto reverse_pass(ColVec&& adj) const {
    ref_type_t<ColVec> adj_ref(std::forward<ColVec>(adj));
    const auto x_size = exp_x_.size();
    Eigen::VectorXd adj_times_jac(x_size);
    double rolling_adjoint_sum = 0.0;
    for (int n = x_size; --n >= 0;) {
      rolling_adjoint_sum += adj_ref(n);
      adj_times_jac(n) = exp_x_(n) * rolling_adjoint_sum;
    }
    return std::make_tuple(adj_times_jac);
  }
};
}  // namespace internal

/**
 * Return an increasing positive ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @tparam ColVec An eigen type with compile time dynamic rows and 1 column
 * @param x Free vector of scalars
 * @return Positive, increasing ordered vector
 */
template <typename ColVec,
          require_eigen_col_vector_vt<is_var, ColVec>* = nullptr>
inline auto positive_ordered_constrain(ColVec&& x) {
  return adj_jac_apply<internal::positive_ordered_constrain_op<ColVec>>(
      std::forward<ColVec>(x));
}

}  // namespace math
}  // namespace stan
#endif
