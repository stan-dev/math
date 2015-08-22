#ifndef STAN_MATH_PRIM_MAT_FUN_INV_PHI_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_PHI_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/inv_Phi.hpp>

namespace stan {
  namespace math {

    /**
     * Return the element-wise standard normal inverse CDF
     *
     * @param p The vector / matrix of cumulative probabilities
     * @return ret(i, j) = inv_Phi(m(i, j))
     */
    template<typename T, int Rows, int Cols>
    inline Eigen::Matrix<T, Rows, Cols>
    inv_Phi(const Eigen::Matrix<T, Rows, Cols>& m) {
      using stan::math::inv_Phi;
      Eigen::Matrix<T, Rows, Cols> ret(m.rows(), m.cols());
      for (size_t j = 0; j < m.cols(); j++)
        for (size_t i = 0; i < m.rows(); i++)
          ret.coeffRef(i,j) = inv_Phi(m.coeffRef(i,j));

      return ret;
    }

  }
}
#endif
