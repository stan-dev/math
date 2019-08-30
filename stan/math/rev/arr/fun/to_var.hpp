#ifndef STAN_MATH_REV_ARR_FUN_TO_VAR_HPP
#define STAN_MATH_REV_ARR_FUN_TO_VAR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/to_var.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] x A std::vector<double>
 * @return A std::vector<var> with the values set
 */
inline std::vector<var> to_var(const std::vector<double>& x) {
  std::vector<var> var_vector(x.size());
  for (size_t n = 0; n < x.size(); n++)
    var_vector[n] = x[n];
  return var_vector;
}

/**
 * Specialization of to_var to pass through pre-constructed vars.
 *
 * Returns the var variable from the input
 *
 * @param[in] x A std::vector<var>
 * @return The input std::vector<var>
 */
template <typename T, enable_if_std_vector<T>...,
          enable_if_var<scalar_type_t<T>>...>
inline auto&& to_var(T&& x) {
  return std::forward<T>(x);
}

}  // namespace math
}  // namespace stan
#endif
