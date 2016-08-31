#ifndef STAN_MATH_PRIM_MAT_FUN_COV_EXP_QUAD_HPP
#define STAN_MATH_PRIM_MAT_FUN_COV_EXP_QUAD_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Returns a squared exponential kernel.
     *
     * @tparam T_x type of std::vector of elements
     * @tparam T_sigma type of sigma
     * @tparam T_l type of length scale
     *
     * @param x std::vector of elements that can be used in square distance.
     *    This function assumes each element of x is the same size.
     * @param sigma standard deviation
     * @param l length scale
     * @return squared distance
     * @throw std::domain_error if sigma <= 0, l <= 0, or
     *   x is nan or infinite
     */
    template<typename T_x, typename T_sigma, typename T_l>
    inline typename
    boost::enable_if_c<!boost::is_same<typename scalar_type<T_x>::type,
                                       double>::value
                       || (boost::is_same<typename scalar_type<T_x>::type,
                                          double>::value
                           & boost::is_same<T_sigma, double>::value
                           & boost::is_same<T_l, double>::value),
    Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                  Eigen::Dynamic, Eigen::Dynamic> >::type
    cov_exp_quad(const std::vector<T_x>& x,
                 T_sigma& sigma,
                 T_l& l) {
      using std::exp;
      check_positive("cov_exp_quad", "sigma", sigma);
      check_positive("cov_exp_quad", "l", l);
      for (size_t n = 0; n < x.size(); n++)
        check_not_nan("cov_exp_quad", "x", x[n]);

      Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                    Eigen::Dynamic, Eigen::Dynamic>
        cov(x.size(), x.size());

      if (x.size() == 0)
        return cov;

      T_sigma sigma_sq = square(sigma);
      T_l neg_half_inv_l_sq = - 0.5 / square(l);

      for (size_t i = 0; i < x.size(); ++i) {
        cov(i, i) = sigma_sq;
        for (size_t j = i + 1; j < x.size(); ++j) {
          cov(i, j) = sigma_sq * exp(squared_distance(x[i], x[j])
                                     * neg_half_inv_l_sq);
          cov(j, i) = cov(i, j);
        }
      }
      return cov;
    }

    /**
     * Returns a squared exponential kernel.
     *
     * @tparam T_x1 type of first std::vector of elements
     * @tparam T_x2 type of second std::vector of elements
     * @tparam T_sigma type of sigma
     * @tparam T_l type of length scale
     *
     * @param x1 std::vector of elements that can be used in square distance
     * @param x2 std::vector of elements that can be used in square distance
     * @param sigma standard deviation
     * @param l length scale
     * @return squared distance
     * @throw std::domain_error if sigma <= 0, l <= 0, or
     *   x is nan or infinite
     */
    template<typename T_x1, typename T_x2, typename T_sigma, typename T_l>
    inline typename
    Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_sigma, T_l>::type,
                  Eigen::Dynamic, Eigen::Dynamic>
    cov_exp_quad(const std::vector<T_x1>& x1,
                 const std::vector<T_x2>& x2,
                 T_sigma& sigma,
                 T_l& l) {
      using std::exp;
      check_positive("cov_exp_quad", "sigma", sigma);
      check_positive("cov_exp_quad", "l", l);
      for (size_t n = 0; n < x1.size(); n++)
        check_not_nan("cov_exp_quad", "x1", x1[n]);
      for (size_t n = 0; n < x2.size(); n++)
        check_not_nan("cov_exp_quad", "x2", x2[n]);

      Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_sigma, T_l>::type,
                    Eigen::Dynamic, Eigen::Dynamic>
        cov(x1.size(), x2.size());
      if (x1.size() == 0 || x2.size() == 0)
        return cov;

      T_sigma sigma_sq = square(sigma);
      T_l neg_half_inv_l_sq = - 0.5 / square(l);

      for (size_t i = 0; i < x1.size(); ++i) {
        for (size_t j = 0; j < x2.size(); ++j) {
          cov(i, j) = sigma_sq * exp(squared_distance(x1[i], x2[j])
                                     * neg_half_inv_l_sq);
        }
      }
      return cov;
    }

  }
}
#endif
