#ifndef STAN_MATH_REV_FUN_INITIALIZE_VARIABLE_HPP
#define STAN_MATH_REV_FUN_INITIALIZE_VARIABLE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Initialize variable to value.  (Function may look pointless, but
 * it's needed to bottom out recursion.)
 */
inline void initialize_variable(var& variable, const var& value) {
  variable = value;
}

/**
 * Initialize every cell in the matrix to the specified value.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 */
template <int R, int C>
inline void initialize_variable(Eigen::Matrix<var, R, C>& matrix,
                                const var& value) {
  matrix.fill(value);
}

/**
 * Initialize the variables in the standard vector recursively.
 *
 * @tparam T type of elements in the vector
 */
template <typename T>
inline void initialize_variable(std::vector<T>& variables, const var& value) {
  for (size_t i = 0; i < variables.size(); ++i) {
    initialize_variable(variables[i], value);
  }
}

}  // namespace math
}  // namespace stan

#endif
