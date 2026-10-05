#ifndef STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP
#define STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_st_var<EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>
eigenvectors_sym(EigMat&& m) {
  if (unlikely(m.size() == 0)) {
    return Eigen::Matrix<value_type_t<EigMat>, -1, -1>(0, 0);
  }
  using PlainMat = plain_type_t<EigMat>;
  decltype(auto) m_ref = to_ref(std::forward<EigMat>(m));
  check_symmetric("eigenvalues_sym", "m", m_ref);

  Eigen::SelfAdjointEigenSolver<PlainMat> solver(
      std::forward<decltype(m_ref)>(m_ref));
  return solver.eigenvectors();
}

}  // namespace math
}  // namespace stan
#endif
