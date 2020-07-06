#ifndef STAN_MATH_REV_CORE_ZERO_ADJOINTS_HPP
#define STAN_MATH_REV_CORE_ZERO_ADJOINTS_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

inline void zero_adjoints();

template <typename T, typename... Pargs, require_st_arithmetic<T>* = nullptr>
inline void zero_adjoints(T& x, Pargs&... args);

template <typename... Pargs>
inline void zero_adjoints(var& x, Pargs&... args);

template <int R, int C, typename... Pargs>
inline void zero_adjoints(Eigen::Matrix<var, R, C>& x, Pargs&... args);

template <typename T, typename... Pargs, require_st_autodiff<T>* = nullptr>
inline void zero_adjoints(std::vector<T>& x, Pargs&... args);

/**
 * End of recursion for set_zero_adjoints
 */
inline void zero_adjoints() {}

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
template <typename T, typename... Pargs, require_st_arithmetic<T>*>
inline void zero_adjoints(T& x, Pargs&... args) {
  zero_adjoints(args...);
}

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
template <typename... Pargs>
inline void zero_adjoints(var& x, Pargs&... args) {
  x.vi_->set_zero_adjoint();
  zero_adjoints(args...);
}

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
template <int R, int C, typename... Pargs>
inline void zero_adjoints(Eigen::Matrix<var, R, C>& x, Pargs&... args) {
  for (size_t i = 0; i < x.size(); ++i)
    x.coeffRef(i).vi_->set_zero_adjoint();
  zero_adjoints(args...);
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
template <typename T, typename... Pargs, require_st_autodiff<T>*>
inline void zero_adjoints(std::vector<T>& x, Pargs&... args) {
  for (size_t i = 0; i < x.size(); ++i)
    zero_adjoints(x[i]);
  zero_adjoints(args...);
}

}  // namespace math
}  // namespace stan
#endif
