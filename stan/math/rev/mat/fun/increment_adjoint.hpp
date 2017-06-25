#ifndef STAN_MATH_REV_MAT_FUN_INCREMENT_ADJOINT_HPP
#define STAN_MATH_REV_MAT_FUN_INCREMENT_ADJOINT_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
  namespace math {

    // A(i, j) must be a var, A.cols(), A.rows() exist
    // B(i, j) must be a double and B must be same size as A
    template <typename T, typename U>
    inline void increment_adjoint(const T& A, const U& B) {
      for (int j = 0; j < A.cols(); ++j)
        for (int i = 0; i < A.rows(); ++i)
          A(i, j).vi_->adj_ += B(i, j);
    }

  }
}
#endif

