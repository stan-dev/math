#ifndef SYMMETRC_EIGENSOLVER_CL_HPP
#define SYMMETRC_EIGENSOLVER_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/mrrr.hpp>
#include <stan/math/opencl/tridiagonalization.hpp>

namespace stan {
namespace math {

template <bool need_eigenvectors = true>
void symmetric_eigensolver(const matrix_cl<double>& A,
                           matrix_cl<double>& eigenvalues,
                           matrix_cl<double>& eigenvectors) {
  matrix_cl<double> packed;
  if (A.rows() <= 2) {
    packed = A;
  } else {
    internal::block_householder_tridiag_cl(A, packed);
  }
  matrix_cl<double> diag = diagonal(packed);
  matrix_cl<double> subdiag = diagonal(
      block_zero_based(packed, 0, 1, packed.rows() - 1, packed.cols() - 1));
  internal::tridiagonal_eigensolver_cl<need_eigenvectors>(
      diag, subdiag, eigenvalues, eigenvectors);
  if (need_eigenvectors && A.rows() > 2) {
    internal::block_apply_packed_Q_cl(packed, eigenvectors);
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
