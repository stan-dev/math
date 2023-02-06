#ifndef STAN_MATH_PRIM_FUN_EIGENVECTORS_HPP
#define STAN_MATH_PRIM_FUN_EIGENVECTORS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Return the eigenvectors of a (real-valued) matrix
 *
 * @tparam EigMat type of real matrix argument
 * @param[in] m matrix to find the eigenvectors of. Must be square and have a
 * non-zero size.
 * @return a complex-valued matrix with columns representing the eigenvectors of
 * `m`
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
          require_not_vt_complex<EigMat>* = nullptr>
inline Eigen::Matrix<complex_return_t<value_type_t<EigMat>>, -1, -1>
eigenvectors(const EigMat& m) {
  using PlainMat = plain_type_t<EigMat>;
  const PlainMat& m_eval = m;
  check_nonzero_size("eigenvectors", "m", m_eval);
  check_square("eigenvectors", "m", m_eval);

  Eigen::EigenSolver<PlainMat> solver(m_eval);
  return solver.eigenvectors();
}

/**
 * Return the eigenvectors of a (complex-valued) matrix
 *
 * @tparam EigCplxMat type of complex matrix argument
 * @param[in] m matrix to find the eigenvectors of. Must be square and have a
 * non-zero size.
 * @return a complex-valued matrix with columns representing the eigenvectors of
 * `m`
 */
template <typename EigCplxMat,
          require_eigen_matrix_dynamic_t<EigCplxMat>* = nullptr,
          require_vt_complex<EigCplxMat>* = nullptr>
inline Eigen::Matrix<complex_return_t<value_type_t<EigCplxMat>>, -1, -1>
eigenvectors(const EigCplxMat& m) {
  using PlainMat = Eigen::Matrix<scalar_type_t<EigCplxMat>, -1, -1>;
  const PlainMat& m_eval = m;
  check_nonzero_size("eigenvectors", "m", m_eval);
  check_square("eigenvectors", "m", m_eval);

  Eigen::ComplexEigenSolver<PlainMat> solver(m_eval);
  return solver.eigenvectors();
}

}  // namespace math
}  // namespace stan
#endif
