#ifndef STAN_MATH_REV_CORE_ZERO_ADJOINTS_HPP
#define STAN_MATH_REV_CORE_ZERO_ADJOINTS_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

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
inline void zero_adjoints(T& x) {
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
inline void zero_adjoints(var& x) {
  x.vi_->adj_ = 0;
}

template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline void zero_adjoints(var_value<EigMat>& x) {
  x.vi_->adj_.setZero();
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
template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
inline void zero_adjoints(EigMat& x) {
  for (size_t i = 0; i < x.size(); ++i) {
    x.coeffRef(i).vi_->adj_ = 0;
  }
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
template <typename T, require_st_autodiff<T>* = nullptr>
inline void zero_adjoints(std::vector<T>& x) {
  for (auto& x_iter : x) {
    zero_adjoints(x_iter);
  }
}

}  // namespace math
}  // namespace stan
#endif
