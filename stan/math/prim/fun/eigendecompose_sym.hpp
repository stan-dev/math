#ifndef STAN_MATH_PRIM_FUN_EIGENDECOMPOSE_SYM_HPP
#define STAN_MATH_PRIM_FUN_EIGENDECOMPOSE_SYM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the eigendecomposition of the specified symmetric matrix.
 *
 * @tparam EigMat type of the matrix
 * @param m Specified matrix.
 * @return A tuple V,D where V is a matrix where the columns are the
 * eigenvectors of m, and D is a column vector of the eigenvalues of m.
 * The eigenvalues are in ascending order of magnitude, with the eigenvectors
 * provided in the same order.
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_st_var<EigMat>* = nullptr>
std::tuple<Eigen::Matrix<value_type_t<EigMat>, -1, -1>,
           Eigen::Matrix<value_type_t<EigMat>, -1, 1>>
eigendecompose_sym(const EigMat& m) {
  if (unlikely(m.size() == 0)) {
    return std::make_tuple(Eigen::Matrix<value_type_t<EigMat>, -1, -1>(0, 0),
                           Eigen::Matrix<value_type_t<EigMat>, -1, 1>(0, 1));
  }
  using PlainMat = plain_type_t<EigMat>;
  const PlainMat& m_eval = m;
  check_symmetric("eigendecompose_sym", "m", m_eval);

  Eigen::SelfAdjointEigenSolver<PlainMat> solver(m_eval);
  return std::make_tuple(std::move(solver.eigenvectors()),
                         std::move(solver.eigenvalues()));
}

}  // namespace math
}  // namespace stan
#endif
