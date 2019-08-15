#ifndef STAN_MATH_REV_MAT_FUN_LOG_DETERMINANT_HPP
#define STAN_MATH_REV_MAT_FUN_LOG_DETERMINANT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

template <int R, int C>
inline var log_determinant(const Eigen::Matrix<var, R, C>& m) {
  using Eigen::Matrix;

  math::check_square("log_determinant", "m", m);

  Eigen::FullPivHouseholderQR<Matrix<double, R, C> > hh
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
