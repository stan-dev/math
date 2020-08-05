#ifndef STAN_MATH_REV_FUN_SOFTMAX_HPP
#define STAN_MATH_REV_FUN_SOFTMAX_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/adj_jac_apply.hpp>
#include <stan/math/rev/functor/ad_stack_matrix.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the softmax of the specified Eigen vector.  Softmax is
 * guaranteed to return a simplex.
 *
 * @param alpha Unconstrained input vector.
 * @return Softmax of the input.
 * @throw std::domain_error If the input vector is size 0.
 */

inline Eigen::Matrix<var, Eigen::Dynamic, 1> softmax(
    const Eigen::Matrix<var, Eigen::Dynamic, 1>& alpha) {
  check_nonzero_size("softmax", "alpha", alpha);

  AD_stack_matrix<Eigen::VectorXd> res_val = softmax(value_of(alpha));
  AD_stack_matrix<Eigen::Matrix<var, Eigen::Dynamic, 1>> res = res_val;
  AD_stack_matrix<Eigen::Matrix<var, Eigen::Dynamic, 1>> alpha_ad_stack = alpha;

  adjac_apply([=]() mutable {
    Eigen::VectorXd res_adj = res.adj();
    alpha_ad_stack.adj() = -res_val * res_adj.dot(res_val) + res_val.cwiseProduct(res_adj);
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
