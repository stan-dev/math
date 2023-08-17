#ifndef STAN_MATH_PRIM_FUN_EIGENDECOMPOSE_HPP
#define STAN_MATH_PRIM_FUN_EIGENDECOMPOSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Return the eigendecomposition of a (real-valued) matrix
 *
 * @tparam EigMat type of real matrix argument
 * @param[in] m matrix to find the eigendecomposition of. Must be square and
 * have a non-zero size.
 * @return A tuple V,D where V is a matrix where the columns are the
 * complex-valued eigenvectors of `m` and D is a complex-valued column vector
 * with entries the eigenvectors of `m`
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
          require_not_vt_complex<EigMat>* = nullptr>
inline std::tuple<Eigen::Matrix<complex_return_t<value_type_t<EigMat>>, -1, -1>,
                  Eigen::Matrix<complex_return_t<value_type_t<EigMat>>, -1, 1>>
eigendecompose(const EigMat& m) {
  if (unlikely(m.size() == 0)) {
    return std::make_tuple(
        Eigen::Matrix<complex_return_t<value_type_t<EigMat>>, -1, -1>(0, 0),
        Eigen::Matrix<complex_return_t<value_type_t<EigMat>>, -1, 1>(0, 1));
  }
  check_square("eigendecompose", "m", m);

  using PlainMat = plain_type_t<EigMat>;
  const PlainMat& m_eval = m;

  Eigen::EigenSolver<PlainMat> solver(m_eval);
  return std::make_tuple(std::move(solver.eigenvectors()),
                         std::move(solver.eigenvalues()));
}

/**
 * Return the eigendecomposition of a (complex-valued) matrix
 *
 * @tparam EigCplxMat type of complex matrix argument
 * @param[in] m matrix to find the eigendecomposition of. Must be square and
 * have a non-zero size.
 * @return A tuple V,D where V is a matrix where the columns are the
 * complex-valued eigenvectors of `m` and D is a complex-valued column vector
 * with entries the eigenvectors of `m`
 */
template <typename EigCplxMat,
          require_eigen_matrix_dynamic_vt<is_complex, EigCplxMat>* = nullptr>
inline std::tuple<
    Eigen::Matrix<complex_return_t<value_type_t<EigCplxMat>>, -1, -1>,
    Eigen::Matrix<complex_return_t<value_type_t<EigCplxMat>>, -1, 1>>
eigendecompose(const EigCplxMat& m) {
  if (unlikely(m.size() == 0)) {
    return std::make_tuple(
        Eigen::Matrix<complex_return_t<value_type_t<EigCplxMat>>, -1, -1>(0, 0),
        Eigen::Matrix<complex_return_t<value_type_t<EigCplxMat>>, -1, 1>(0, 1));
  }
  check_square("eigendecompose", "m", m);

  using PlainMat = Eigen::Matrix<scalar_type_t<EigCplxMat>, -1, -1>;
  const PlainMat& m_eval = m;

  Eigen::ComplexEigenSolver<PlainMat> solver(m_eval);

  return std::make_tuple(std::move(solver.eigenvectors()),
                         std::move(solver.eigenvalues()));
}

}  // namespace math
}  // namespace stan
#endif
