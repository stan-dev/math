#ifndef STAN_MATH_PRIM_FUN_EIGENVALUES_HPP
#define STAN_MATH_PRIM_FUN_EIGENVALUES_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

template <typename T, require_eigen_t<T>* = nullptr>
inline auto eigenvalues(const T& m) {
  check_nonzero_size("eigenvalues", "m", m);
  check_square("eigenvalues", "m", m);
  return Eigen::EigenSolver<plain_type_t<T>>(m).eigenvalues().eval();
}

}  // namespace math
}  // namespace stan
#endif
