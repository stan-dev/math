#ifndef STAN_MATH_PRIM_FUN_EVAL_HPP
#define STAN_MATH_PRIM_FUN_EVAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Inputs that are not Eigen types are forwarded unmodified
 *
 * @tparam T Input type
 * @param[in] arg Input argument
 * @return Forwarded input argument
 **/
template <typename T,
	  require_not_eigen_t<T>* = nullptr>
inline decltype(auto) eval(T&& arg) {
  return std::forward<T>(arg);
}

/**
 * Inputs that are Eigen types are eval'd and returned
 *
 * @tparam T Input type
 * @param[in] arg Input argument
 * @return Forwarded input argument
 **/
template <typename T,
	  require_eigen_t<T>* = nullptr>
inline decltype(auto) eval(T&& arg) {
  return arg.eval();
}

}  // namespace math
}  // namespace stan

#endif
