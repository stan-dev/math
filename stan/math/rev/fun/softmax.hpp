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
class softmax_op {
  int N_;
  double* y_;  // Holds the results of the softmax

 public:
  softmax_op() : N_(0), y_(nullptr) {}

  /**
   * Compute the softmax of the unconstrained input vector
   *
   * @param alpha Unconstrained input vector.
   * @return Softmax of the input.
   */
  template <typename T, std::size_t size>
  plain_type_t<T> operator()(const std::array<bool, size>& /* needs_adj */,
                     const T& alpha) {
    N_ = alpha.size();
    y_ = ChainableStack::instance_->memalloc_.alloc_array<double>(N_);
    Eigen::Map<plain_type_t<T>> y_map(y_, N_);

    y_map = softmax(alpha);

    return y_map;
  }

  /**
   * Compute the result of multiply the transpose of the adjoint vector times
   * the Jacobian of the softmax operator. It is more efficient to do this
   * without actually computing the Jacobian and doing the vector-matrix
   * product.
   *
   * @param adj Eigen::VectorXd of adjoints at the output of the softmax
   * @return Eigen::VectorXd of adjoints propagated through softmax operation
   */
  template <typename T, std::size_t size, typename T_plain = plain_type_t<T>>
  std::tuple<T_plain> multiply_adjoint_jacobian(
      const std::array<bool, size>& /* needs_adj */,
      const T& adj) const {
    T_plain adj_times_jac(N_);
    Eigen::Map<T_plain> y(y_, N_);

    adj_times_jac = -y * adj.dot(y) + y.cwiseProduct(adj);

    return std::make_tuple(adj_times_jac);
  }
};
}  // namespace internal

/**
 * Return the softmax of the specified Eigen vector.  Softmax is
 * guaranteed to return a simplex.
 *
 * @param alpha Unconstrained input vector.
 * @return Softmax of the input.
 * @throw std::domain_error If the input vector is size 0.
 */
template <typename T, require_t<is_var<scalar_type_t<T>>>...>
inline auto softmax(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& alpha) {
    const Eigen::Ref<const plain_type_t<T>>& a_ref = alpha;
    check_nonzero_size("softmax", "alpha", alpha);

    return adj_jac_apply<internal::softmax_op>(a_ref);
  });
}

}  // namespace math
}  // namespace stan
#endif
