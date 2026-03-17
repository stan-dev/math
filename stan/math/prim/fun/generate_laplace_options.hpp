#ifndef STAN_MATH_PRIM_FUN_GENERATE_LAPLACE_OPTIONS_HPP
#define STAN_MATH_PRIM_FUN_GENERATE_LAPLACE_OPTIONS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

namespace internal {
inline constexpr int laplace_default_hessian_block_size = 1;
inline constexpr int laplace_default_solver = 1;
inline constexpr double laplace_default_tolerance = 1.49012e-08;
inline constexpr int laplace_default_max_num_steps = 500;
inline constexpr int laplace_default_allow_fallthrough = 1;
inline constexpr int laplace_default_max_steps_line_search = 1000;
}  // namespace internal

/**
 * User function for generating laplace options tuple
 * @param theta_0_size Size of user supplied initial theta
 * @return tuple representing laplace options exposed to user.
 */
inline auto generate_laplace_options(int theta_0_size) {
  return std::make_tuple(Eigen::VectorXd::Zero(theta_0_size).eval(),
                         internal::laplace_default_tolerance,
                         internal::laplace_default_max_num_steps,
                         internal::laplace_default_solver,
                         internal::laplace_default_max_steps_line_search,
                         internal::laplace_default_allow_fallthrough);
}

/**
 * User function for generating laplace options tuple
 * @tparam ThetaVec An Eigen vector type for user supplied initial theta
 * @param theta_0 User supplied initial theta
 * @return tuple representing laplace options exposed to user.
 */
template <typename ThetaVec, require_eigen_t<ThetaVec>* = nullptr>
inline auto generate_laplace_options(ThetaVec&& theta_0) {
  return std::make_tuple(stan::math::eval(std::forward<ThetaVec>(theta_0)),
                         internal::laplace_default_tolerance,
                         internal::laplace_default_max_num_steps,
                         internal::laplace_default_solver,
                         internal::laplace_default_max_steps_line_search,
                         internal::laplace_default_allow_fallthrough);
}

}  // namespace math
}  // namespace stan

#endif
