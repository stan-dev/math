#ifndef STAN_MATH_PRIM_MAT_FUN_COV_SQ_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_COV_SQ_EXP_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/square.hpp>

namespace stan {
  namespace math {
    /**
     * Returns a squared exponential kernel.
     *
     * @param x Vector or Matrix
     * @param sigma standard deviation
     * @param l length scale
     * @return Dot product of the vectors.
     * @throw std::domain_error If the vectors are not the same
     * size or if they are both not vector dimensioned.
     */
    template<int R1, int C1, typename T_x, typename T_sigma, typename T_l>
    inline typename
    Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                  Eigen::Dynamic, Eigen::Dynamic>
    cov_sq_exp(const Eigen::Matrix<T_x, R1, C1>& x,
               T_sigma& sigma,
               T_l& l) {
      stan::math::check_positive("cov_sq_exp", "sigma", sigma);
      stan::math::check_positive("cov_sq_exp", "l", l);
      stan::math::check_not_nan("cov_sq_exp", "x", x);

      Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                    Eigen::Dynamic, Eigen::Dynamic>
        cov(x.rows(), x.rows());

      T_sigma sigma_sq = square(sigma);
      T_l neg_half_inv_l_sq = - 0.5 / square(l);

      
      
      for (int i = 0; i < x.rows(); i++) {
        for (int j = i + 1; j < x.rows(); j++) {
          // sigma 
          cov(i, j) = sigma_sq * exp((x.row(i) - x.row(j)).squaredNorm() * neg_half_inv_l_sq);
          cov(j, i) = cov(i, j);
        }
      }
      for (int i = 0; i < x.rows(); i++) {
        cov(i, i) = sigma_sq;
      }
      return cov;
    }
  }
}
#endif
