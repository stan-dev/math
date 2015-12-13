#ifndef STAN_MATH_PRIM_MAT_FUN_COV_SQ_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_COV_SQ_EXP_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/square.hpp>

namespace stan {
  namespace math {

    /**
     * Returns a squared exponential kernel.
     *
     * @param x std::vector of a type that can be square distance
     * @param sigma standard deviation
     * @param l length scale
     * @return squared distance
     * @throw std::domain_error if sigma <= 0, l <= 0, or 
     *   x is nan or infinite
     */
    template<typename T_x, typename T_sigma, typename T_l>
    inline typename
    Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                  Eigen::Dynamic, Eigen::Dynamic>
    cov_sq_exp(const std::vector<T_x>& x,
               T_sigma& sigma,
               T_l& l) {
      using std::exp;
      stan::math::check_positive("cov_sq_exp", "sigma", sigma);
      stan::math::check_positive("cov_sq_exp", "l", l);
      stan::math::check_not_nan("cov_sq_exp", "x", x);

      Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                    Eigen::Dynamic, Eigen::Dynamic>
        cov(x.size(), x.size());

      T_sigma sigma_sq = square(sigma);
      T_l neg_half_inv_l_sq = - 0.5 / square(l);
      
      for (size_t i = 0; i < x.size(); i++) {
        for (size_t j = i + 1; j < x.size(); j++) {
          cov(i, j) = sigma_sq * exp(squared_distance(x[i], x[j]) * neg_half_inv_l_sq);
          cov(j, i) = cov(i, j);
        }
      }
      for (size_t i = 0; i < x.size(); i++)
        cov(i, i) = sigma_sq;
      return cov;
    }
  }
}
#endif
