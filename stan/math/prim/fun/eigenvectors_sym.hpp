#ifndef STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP
#define STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP

#include <stanh/prim/err/check_nonzero_size.hpp>
#include <stanh/prim/err/check_symmetric.hpp>
#include <stanh/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigenvectors_sym(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
  check_nonzero_size("eigenvectors_sym", "m", m);
  check_symmetric("eigenvalues_sym", "m", m);

  Eigen::SelfAdjointEigenSolver<
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
      solver(m);
  return solver.eigenvectors();
}

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP
#define STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP

#include <stanh/prim/err/check_nonzero_size.hpp>
#include <stanh/prim/err/check_symmetric.hpp>
#include <stanh/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigenvectors_sym(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
  check_nonzero_size("eigenvectors_sym", "m", m);
  check_symmetric("eigenvalues_sym", "m", m);

  Eigen::SelfAdjointEigenSolver<
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
      solver(m);
  return solver.eigenvectors();
}

}  // namespace math
}  // namespace stan
#endif
