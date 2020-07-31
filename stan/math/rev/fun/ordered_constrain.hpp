#ifndef STAN_MATH_REV_FUN_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_ORDERED_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/adj_jac_apply.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename T>
class ordered_constrain_op {
 public:
  adj_arg<T> exp_x_;
  const size_t x_size_;
  explicit ordered_constrain_op(const T& x)
      : exp_x_(x.size() > 0 ? x.size() - 1 : 0), x_size_(x.size()) {}
  /**
   * Return an increasing ordered vector derived from the specified
   * free vector.  The returned constrained vector will have the
   * same dimensionality as the specified free vector.
   *
   * @tparam ColVec An eigen type with compile time dynamic rows and 1 column
   * @param x Free vector of scalars
   * @return Increasing ordered vector
   */
  template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
  Eigen::VectorXd operator()(ColVec&& x) {
    using std::exp;
    ref_type_t<ColVec> x_ref(std::forward<ColVec>(x));
    Eigen::Matrix<double, Eigen::Dynamic, 1> y(x_size_);
    if (x_size_ == 0) {
      return y;
    }
    y[0] = x_ref[0];
    for (int n = 1; n < x_size_; ++n) {
      exp_x_(n - 1) = exp(x_ref[n]);
      y[n] = y[n - 1] + exp_x_(n - 1);
    }
    return y;
  }

  /**
   * Compute the result of multiply the transpose of the adjoint vector times
   * the Jacobian of the ordered_constrain operator.
   *
   * @tparam ColVec An eigen type with compile time dynamic rows and 1 column
   * @param adj Eigen::VectorXd of adjoints at the output of the softmax
   * @return Eigen::VectorXd of adjoints propagated through softmax operation
   */
  template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
  auto multiply_adjoint_jacobian(ColVec&& adj) const {
  ref_type_t<ColVec> adj_ref(std::forward<ColVec>(adj));
  Eigen::VectorXd adj_times_jac(x_size_);
    double rolling_adjoint_sum = 0.0;
    if (x_size_ > 0) {
      for (int n = x_size_ - 1; n > 0; --n) {
        rolling_adjoint_sum += adj_ref(n);
        adj_times_jac(n) = exp_x_(n - 1) * rolling_adjoint_sum;
      }
      adj_times_jac(0) = rolling_adjoint_sum + adj_ref(0);
    }
    return std::make_tuple(adj_times_jac);
  }
};
}  // namespace internal

/**
 * Return an increasing ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @tparam ColVec An eigen type with compile time dynamic rows and 1 column
 * @param x Free vector of scalars
 * @return Increasing ordered vector
 */
template <typename ColVec, require_eigen_col_vector_vt<is_var, ColVec>* = nullptr>
inline auto ordered_constrain(ColVec&& x) {
  return adj_jac_apply<internal::ordered_constrain_op<ColVec>>(std::forward<ColVec>(x));
}

}  // namespace math
}  // namespace stan
#endif
