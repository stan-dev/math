#ifndef STAN_MATH_REV_CORE_ZERO_ADJOINTS_HPP
#define STAN_MATH_REV_CORE_ZERO_ADJOINTS_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * End of recursion for set_zero_adjoints
 */
inline void zero_adjoints() noexcept {}

/**
 * Do nothing for non-autodiff arguments. Recursively call zero_adjoints
 * on the rest of the arguments.
 *
 * @tparam T type of current argument
 * @tparam Pargs type of rest of arguments
 *
 * @param x current argument
 * @param args rest of arguments to zero
 */
template <typename T, require_st_arithmetic<T>* = nullptr>
inline void zero_adjoints(T& x) noexcept {}

/**
 * Zero the adjoint of the vari in the first argument. Recursively call
 * zero_adjoints on the rest of the arguments.
 *
 * @tparam T type of current argument
 * @tparam Pargs type of rest of arguments
 *
 * @param x current argument
 * @param args rest of arguments to zero
 */
inline void zero_adjoints(var& x) { x.adj() = 0; }

/**
 * Zero the adjoints of the varis of every var in an Eigen::Matrix
 * container. Recursively call zero_adjoints on the rest of the arguments.
 *
 * @tparam T type of current argument
 * @tparam Pargs type of rest of arguments
 *
 * @param x current argument
 * @param args rest of arguments to zero
 */
template <typename EigMat, require_eigen_vt<is_autodiff, EigMat>* = nullptr>
inline void zero_adjoints(EigMat& x) {
  for (size_t i = 0; i < x.size(); ++i)
    x.coeffRef(i).adj() = 0;
}

/**
 * Zero the adjoints of every element in a vector. Recursively call
 * zero_adjoints on the rest of the arguments.
 *
 * @tparam T type of current argument
 * @tparam Pargs type of rest of arguments
 *
 * @param x current argument
 * @param args rest of arguments to zero
 */
template <typename StdVec,
          require_std_vector_st<is_autodiff, StdVec>* = nullptr>
inline void zero_adjoints(StdVec& x) {
  for (size_t i = 0; i < x.size(); ++i) {
    zero_adjoints(x[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif
