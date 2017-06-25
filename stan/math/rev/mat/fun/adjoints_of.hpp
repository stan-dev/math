#ifndef STAN_MATH_REV_MAT_FUN_ADJOINTS_OF_HPP
#define STAN_MATH_REV_MAT_FUN_ADJOINTS_OF_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
  namespace math {

    // S must be matrix-like with R row type and C col type
    template <int R, int C, typename S>
    inline Eigen::Matrix<double, R, C> adjoints_of(const S& A) {
      Eigen::Matrix<double, R, C> adjA(A.rows(), A.cols());
      for (int j = 0; j < A.cols(); ++j)
        for (int i = 0; i < A.rows(); ++i)
          adjA(i, j) = A(i, j).vi_->adj_;
      return adjA;
    }

  }
}
#endif
