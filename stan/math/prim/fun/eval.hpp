#ifndef STAN_MATH_PRIM_FUN_EVAL_HPP
#define STAN_MATH_PRIM_FUN_EVAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Inputs which have a plain_type equal to the own time are forwarded
 * unmodified (for Eigen expressions these types are different)
 *
 * @tparam T Input type
 * @param[in] arg Input argument
 * @return Forwarded input argument
 **/
template <typename T,
          require_same_t<std::decay_t<T>, plain_type_t<T>>* = nullptr>
inline T eval(T&& arg) {
  return std::forward<T>(arg);
}

/**
 * Inputs which have a plain_type different from their own type are
 * Eval'd (this catches Eigen expressions)
 *
 * @tparam T Input type
 * @param[in] arg Input argument
 * @return Eval'd argument
 **/
template <typename T,
          require_not_same_t<std::decay_t<T>, plain_type_t<T>>* = nullptr>
inline decltype(auto) eval(const T& arg) {
  return arg.eval();
}

}  // namespace math
}  // namespace stan

#endif
