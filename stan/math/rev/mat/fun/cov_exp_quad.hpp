#ifndef STAN_MATH_REV_MAT_FUN_COV_EXP_QUAD_HPP
#define STAN_MATH_REV_MAT_FUN_COV_EXP_QUAD_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/meta/length.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <vector>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * This is a subclass of the vari class for precomputed
     * gradients of cov_exp_quad.
     *
     * The class stores the double values for the distance
     * matrix, pointers to the varis for the covariance
     * matrix, along with a pointer to the vari for sigma,
     * and the vari for l.
     *
     * @tparam T_x type of std::vector of elements
     * @tparam T_sigma type of sigma
     * @tparam T_l type of length scale
     */
    template <typename T_x, typename T_sigma, typename T_l>
    class cov_exp_quad_vari : public vari {
      public:
        int size_;
        int size_sq_;
        int size_vech_;
        double l_d_;
        double sigma_d_;
        double sigma_sq_d_;
        double* dist_;
        vari* l_vari_;
        vari* sigma_vari_;
        vari** cov_lower_;
        vari** cov_upper_;
        vari** cov_diag_;

        /**
         * Constructor for cov_exp_quad.
         *
         * All memory allocated in
         * ChainableStack's stack_alloc arena.
         *
         * It is critical for the efficiency of this object
         * that the constructor create new varis that aren't
         * popped onto the var_stack_, but rather are
         * popped onto the var_nochain_stack_. This is
         * controlled to the second argument to
         * vari's constructor.
         *
         * @param x std::vector input that can be used in square distance
         *    Assumes each element of x is the same size
         * @param sigma standard deviation
         * @param l length scale
         */
        cov_exp_quad_vari(const std::vector<T_x>& x,
                          const T_sigma& sigma,
                          const T_l& l)
          : vari(0.0),
          size_(x.size()),
          size_vech_(size_ * (size_ - 1) / 2),
          l_d_(value_of(l)), sigma_d_(value_of(sigma)),
          sigma_sq_d_(std::pow(sigma_d_, 2)),
          dist_(ChainableStack::memalloc_.alloc_array<double>(size_vech_)),
          l_vari_(l.vi_), sigma_vari_(sigma.vi_),
          cov_lower_(ChainableStack::memalloc_.alloc_array<vari*>(size_vech_)), 
          cov_upper_(ChainableStack::memalloc_.alloc_array<vari*>(size_vech_)),
          cov_diag_(ChainableStack::memalloc_.alloc_array<vari*>(size_)) {
            double inv_half_sq_l_d = 0.5 / (std::pow(l_d_, 2));
            size_t pos = 0;
            for (size_t j = 0; j < static_cast<size_t>(size_ - 1); ++j)
              for (size_t i = j + 1; i < static_cast<size_t>(size_); ++i) {
                double dist_sq = squared_distance(x[i],x[j]);
                double val = sigma_sq_d_ * exp(-dist_sq
                                              * inv_half_sq_l_d);
                dist_[pos] = dist_sq;
                cov_upper_[pos] = new vari(val, false);
                cov_lower_[pos] = new vari(val, false);
                ++pos;
              }
              for (size_t i = 0; i < static_cast<size_t>(size_); ++i) 
                cov_diag_[i] = new vari(sigma_sq_d_, false);
          }

        virtual void chain() {
          double adjl(0.0);
          double adjsigma(0.0);

          for (size_t i = 0; i < static_cast<size_t>(size_vech_); ++i) {
            vari* el_low = cov_lower_[i];
            vari* el_high = cov_upper_[i];
            double prod_adj_cov_low = el_low->adj_ * el_low->val_;
            double prod_adj_cov_high = el_high->adj_ * el_high->val_;
            double prod_add = prod_adj_cov_low + prod_adj_cov_high;
            adjl += prod_add * dist_[i];
            adjsigma += prod_add;
          }
          for (size_t i = 0; i < static_cast<size_t>(size_); ++i) {
            vari* el = cov_diag_[i];
            double prod_adj_cov = el->adj_ * el->val_;
            adjsigma += prod_adj_cov;
          }
          l_vari_->adj_ +=  adjl / std::pow(l_d_, 3);
          sigma_vari_->adj_ += adjsigma * 2 / sigma_d_;
          }
    };

    /**
     * Returns a squared exponential kernel.
     *
     * @param x std::vector input that can be used in square distance
     *    Assumes each element of x is the same size
     * @param sigma standard deviation
     * @param l length scale
     * @return squared distance
     * @throw std::domain_error if sigma <= 0, l <= 0, or
     *   x is nan or infinite
     */
    template <>
    inline
    Eigen::Matrix<var, -1, -1>
    cov_exp_quad(const std::vector<double>& x,
                 var& sigma,
                 var& l) {
      check_positive("cov_exp_quad", "sigma", sigma);
      check_positive("cov_exp_quad", "l", l);
      for (size_t n = 0; n < x.size(); n++)
        check_not_nan("cov_exp_quad", "x", x[n]);

      int x_size = x.size();
      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
        cov(x_size, x_size);
      if (x_size == 0)
        return cov;

      cov_exp_quad_vari<double, var, var> *baseVari
        = new cov_exp_quad_vari<double, var, var>(x, sigma, l);

      size_t pos = 0;
      for (size_t j = 0; j < static_cast<size_t>(x_size - 1); ++j) {
        for (size_t i = (j + 1); i < static_cast<size_t>(x_size); ++i) {
          cov.coeffRef(i, j).vi_ = baseVari->cov_lower_[pos];
          cov.coeffRef(j, i).vi_ = baseVari->cov_upper_[pos];
          ++pos;
        }
        cov.coeffRef(j, j).vi_ = baseVari->cov_diag_[j];
      }
      cov.coeffRef(x_size - 1, x_size - 1).vi_ = baseVari->cov_diag_[x_size - 1];
      return cov;
    }
  }
}
#endif
