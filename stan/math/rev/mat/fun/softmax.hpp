#ifndef STAN_MATH_REV_MAT_FUN_SOFTMAX_HPP
#define STAN_MATH_REV_MAT_FUN_SOFTMAX_HPP

#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/softmax.hpp>
#include <stan/math/rev/mat.hpp>
#include <vector>

namespace stan {
namespace math {

namespace {
class softmax_op {
  int N;
  double* y_mem_;  // Holds the results of the softmax
 public:
  /*
   * Compute the softmax of the unconstrained input vector
   *
   * @param alpha Unconstrained input vector.
   * @return Softmax of the input.
   */
  Eigen::VectorXd apply(const Eigen::VectorXd& alpha) {
    N = alpha.size();
    y_mem_ = ChainableStack::instance().memalloc_.alloc_array<double>(N);

    auto y = softmax(alpha);
    for (int n = 0; n < N; ++n)
      y_mem_[n] = y(n);
    return y;
  }

  /*
   * Compute the result of multiply the transpose of the adjoint vector times
   * the Jacobian of the softmax operator. It is more efficient to do this
   * without actually computing the Jacobian and doing the vector-matrix
   * product.
   *
   * @param adj Eigen::VectorXd of adjoints at the output of the softmax
   * @return Eigen::VectorXd of adjoints propagated through softmax operation
   */
  Eigen::VectorXd multiply_adjoint_jacobian(const Eigen::VectorXd& adj) const {
    Eigen::VectorXd adj_times_jac(N);

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > y(y_mem_, N);

    double adj_dot_y = adj.dot(y);
    for (int n = 0; n < N; ++n) {
      adj_times_jac(n) = -y(n) * adj_dot_y + y(n) * adj(n);
    }

    return adj_times_jac;
  }

  //  Allocate on the autodiff stack
  static inline void* operator new(size_t nbytes) {
    return ChainableStack::instance().memalloc_.alloc(nbytes);
  }

  //  No need to call the destructor, that'll get handled elsewhere
  static inline void operator delete(void*) {}
};
}  // namespace

/**
 * Return the softmax of the specified Eigen vector.  Softmax is
 * guaranteed to return a simplex.
 *
 * @param alpha Unconstrained input vector.
 * @return Softmax of the input.
 * @throw std::domain_error If the input vector is size 0.
 */
Eigen::Matrix<var, Eigen::Dynamic, 1> softmax(
    const Eigen::Matrix<var, Eigen::Dynamic, 1>& alpha) {
  check_nonzero_size("softmax", "alpha", alpha);

  return adj_jac_apply<softmax_op>(alpha);
}

}  // namespace math
}  // namespace stan
#endif
