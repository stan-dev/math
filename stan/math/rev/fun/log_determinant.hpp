#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
inline var log_determinant(const EigMat& m) {
  using Eigen::Matrix;

  math::check_square("log_determinant", "m", m);

  Eigen::FullPivHouseholderQR<promote_scalar_t<double, EigMat>> hh
      = m.val().fullPivHouseholderQr();

  vari** varis
      = ChainableStack::instance_->memalloc_.alloc_array<vari*>(m.size());
  Eigen::Map<matrix_vi>(varis, m.rows(), m.cols()) = m.vi();

  double* gradients
      = ChainableStack::instance_->memalloc_.alloc_array<double>(m.size());
  Eigen::Map<matrix_d>(gradients, m.rows(), m.cols())
      = hh.inverse().transpose();

  return var(new precomputed_gradients_vari(hh.logAbsDeterminant(), m.size(),
                                            varis, gradients));
}

}  // namespace math
}  // namespace stan
#endif
