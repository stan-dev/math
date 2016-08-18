#ifndef STAN_MATH_REV_MAT_FUN_MATRIX_EXP_HPP
#define STAN_MATH_REV_MAT_FUN_MATRIX_EXP_HPP

#include <unsupported/Eigen/MatrixFunctions>

namespace stan {
  namespace math {

      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
      matrix_exp(const Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>& A){
          
          return A.exp();
      }

  }
}
#endif
