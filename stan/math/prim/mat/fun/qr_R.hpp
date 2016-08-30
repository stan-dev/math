#ifndef STAN_MATH_PRIM_MAT_FUN_QR_R_HPP
#define STAN_MATH_PRIM_MAT_FUN_QR_R_HPP

#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <Eigen/QR>

namespace stan {
  namespace math {

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    qr_R(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
      typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
      check_nonzero_size("qr_R", "m", m);
      check_greater_or_equal("qr_R",
                             "m.rows()",
                             static_cast<size_t>(m.rows()),
                             static_cast<size_t>(m.cols()));
      Eigen::HouseholderQR<matrix_t> qr(m.rows(), m.cols());
      qr.compute(m);
      matrix_t R = qr.matrixQR();
      if (m.rows() > m.cols())
        R.bottomRows(m.rows() - m.cols()).setZero();
      for (int i = 0; i < R.cols(); i++) {
        for (int j = 0; j < i; j++)
          R.coeffRef(i, j) = 0.0;
        if (R(i, i) < 0)
          R.row(i) *= -1.0;
      }
      return R;
    }
  }
}
#endif
