#ifndef STAN_MATH_REV_FUN_SOFTMAX_HPP
#define STAN_MATH_REV_FUN_SOFTMAX_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/adj_jac_apply.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {

template <typename T>
class softmax_op {
  adj_arg_t<T> y_;

 public:
  explicit softmax_op(const T& x) : y_(setup_adj_arg<T>(x.size())) {}

  /**
   * Compute the softmax of the unconstrained input vector
   *
   * @tparam ColVec An eigen type with dynamic rows and 1 column at compile
   * time.
   * @param alpha Unconstrained input vector.
   * @return Softmax of the input.
   */
  template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
  auto operator()(ColVec&& alpha) {
    y_ = softmax(std::forward<ColVec>(alpha));
    return y_.eval();
  }

  /**
   * Compute the result of multiply the transpose of the adjoint vector times
   * the Jacobian of the softmax operator. It is more efficient to do this
   * without actually computing the Jacobian and doing the vector-matrix
   * product.
   *
   * @tparam ColVec An eigen type with dynamic rows and 1 column at compile
   * time.
   * @param adj Eigen::VectorXd of adjoints at the output of the softmax
   * @return Eigen::VectorXd of adjoints propagated through softmax operation
   */
  template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
  auto multiply_adjoint_jacobian(ColVec&& adj) const {
    ref_type_t<ColVec> adj_ref(std::forward<ColVec>(adj));
    return std::make_tuple(-y_ * adj_ref.dot(y_) + y_.cwiseProduct(adj_ref));
  }
};
}  // namespace internal

/**
 * Return the softmax of the specified Eigen vector.  Softmax is
 * guaranteed to return a simplex.
 *
 * @tparam ColVec An eigen type with dynamic rows and 1 column at compile time.
 * @param alpha Unconstrained input vector.
 * @return Softmax of the input.
 * @throw std::domain_error If the input vector is size 0.
 */
template <typename ColVec,
          require_eigen_col_vector_vt<is_var, ColVec>* = nullptr>
inline auto softmax(ColVec&& alpha) {
  check_nonzero_size("softmax", "alpha", alpha);
  return adj_jac_apply<internal::softmax_op<ColVec>>(
      std::forward<ColVec>(alpha));
}

}  // namespace math
}  // namespace stan
#endif
