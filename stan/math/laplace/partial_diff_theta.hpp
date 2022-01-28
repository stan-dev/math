#ifndef STAN_MATH_LAPLACE_PARTIAL_DIFF_THETA_HPP
#define STAN_MATH_LAPLACE_PARTIAL_DIFF_THETA_HPP

// TODO: refine include.
#include <stan/math/mix.hpp>
#include <stan/math/prim/fun/col.hpp>

namespace stan {
namespace math {
/**
 * Returns the partial derivative of the approximate marginal
 * distribution with respect to theta and eta.
 * The derivative with respect to theta is denoted s2 in
 * laplace_marginal.hpp.
 */
// TODO: rename function, since we also differentiate wrt eta.
// TODO: address case where eta / theta are doubles and we don't
// want full derivatives.
template <typename F>
inline Eigen::VectorXd partial_diff_theta(
    const F& f_, const Eigen::VectorXd& theta, const Eigen::VectorXd& eta,
    const Eigen::VectorXd& delta, const std::vector<int>& delta_int,
    const Eigen::MatrixXd& A, int hessian_block_size,
    std::ostream* pstream = 0) {
}

}  // namespace math
}  // namespace stan

#endif
