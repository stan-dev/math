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
#include <stan/math/rev/mat/type/symmetric.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <vector>
#include <iostream>
#include <cmath>

namespace stan {
  namespace math {

    /**
     */
    class cov_exp_quad_vari : public vari {
      public:
        int size_;
        int tri_size_;
        double l_d_;
        double sigma_d_;
        double sigma_sq_d_;
        double* dist_;
        vari* l_vari_;
        vari* sigma_vari_;
        vari** cov_;

        /** 
         * Constructor for multiply_mat_vari.
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
         * @param A matrix
         * @param B matrix
         */
        cov_exp_quad_vari(const std::vector<double>& x,
                          const var& sigma,
                          const var& l)
          : vari(0.0),
          size_(x.size()),
          tri_size_(size_ * (size_ + 1) / 2),
          l_d_(value_of(l)), sigma_d_(value_of(sigma)),
          sigma_sq_d_(std::pow(sigma_d_,2)),
          dist_(ChainableStack::memalloc_.alloc_array<double>(tri_size_)),
          l_vari_(l.vi_), sigma_vari_(sigma.vi_),
          cov_(ChainableStack::memalloc_.alloc_array<vari*>(tri_size_)) {
            size_t pos = 0;
            double inv_half_sq_l_d = 0.5 / (std::pow(l_d_,2));
            for (size_t j = 0; j < static_cast<size_t>(size_); ++j) 
              for (size_t i = j; i < static_cast<size_t>(size_); ++i) {
                dist_[pos] = squared_distance(x[i],x[j]) 
                             * inv_half_sq_l_d;
                cov_[pos] = new vari((i == j) ? sigma_sq_d_ : sigma_sq_d_ * exp(-dist_[pos]), false);
                ++pos;
              }
          }

        virtual void chain() {
          using Eigen::MatrixXd;
          using Eigen::ArrayXd;
          using Eigen::Map;
          double adjl;
          double adjsigma;
          ArrayXd adj_cov(tri_size_);
          ArrayXd cov(tri_size_);

          for (size_t i = 0; i < static_cast<size_t>(tri_size_); ++i) {
            adj_cov(i) = cov_[i]->adj_;
            cov(i) = cov_[i]->val_;
          }
          adjl = (adj_cov * Map<ArrayXd>(dist_, tri_size_)
            * cov).sum() * 2 / l_d_;
          adjsigma = 2 / sigma_d_ * (adj_cov * cov).sum();
          l_vari_->adj_ += adjl;
          sigma_vari_->adj_ += adjsigma;
        }
    };

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
    inline symmetric 
    cov_exp_quad(const std::vector<double>& x,
                 const var& sigma,
                 const var& l) {
      check_positive("cov_exp_quad", "sigma", sigma);
      check_positive("cov_exp_quad", "l", l);
      for (size_t n = 0; n < x.size(); n++)
        check_not_nan("cov_exp_quad", "x", x[n]);

      symmetric cov(x.size());
      if (x.size() == 0)
        return cov;

      cov_exp_quad_vari *baseVari
        = new cov_exp_quad_vari(x, sigma, l);

      size_t pos = 0;
      for (size_t j = 0; j < x.size(); ++j)
        for (size_t i = j; i < x.size(); ++i)
          cov(i,j).vi_ = baseVari->cov_[pos++];
      return cov;
    }
  }
}
#endif
