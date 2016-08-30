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
        double l_d_;
        double sigma_d_;
        double sigma_sq_d_;
        double* dist_;
        vari* l_vari_;
        vari* sigma_vari_;
        vari** cov_;

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
          l_d_(value_of(l)), sigma_d_(value_of(sigma)),
          sigma_sq_d_(std::pow(sigma_d_, 2)),
          dist_(ChainableStack::memalloc_.alloc_array<double>(size_
                                                                  * size_)),
          l_vari_(l.vi_), sigma_vari_(sigma.vi_),
          cov_(ChainableStack::memalloc_.alloc_array<vari*>(size_
                                                            * size_)) {
            size_t pos = 0;
            double inv_half_sq_l_d = 0.5 / (std::pow(l_d_, 2));
            for (size_t j = 0; j < static_cast<size_t>(size_); ++j)
              for (size_t i = 0; i < static_cast<size_t>(size_); ++i) {
                dist_[pos] = squared_distance(x[i], x[j])
                             * inv_half_sq_l_d;
                cov_[pos] = new vari((i == j) ? sigma_sq_d_
                                     : sigma_sq_d_ * exp(-dist_[pos]), false);
                ++pos;
              }
          }

        virtual void chain() {
          using Eigen::MatrixXd;
          using Eigen::ArrayXXd;
          using Eigen::Map;
          double adjl;
          double adjsigma;
          ArrayXXd adj_cov(size_, size_);
          ArrayXXd cov(size_, size_);

          for (size_t i = 0; i < static_cast<size_t>(adj_cov.size()); ++i) {
            adj_cov(i) = cov_[i]->adj_;
            cov(i) = cov_[i]->val_;
          }
          adjl = (adj_cov * Map<ArrayXXd>(dist_, size_, size_)
            * cov).sum() * 2 / l_d_;
          adjsigma = 2 / sigma_d_ * (adj_cov * cov).sum();
          l_vari_->adj_ += adjl;
          sigma_vari_->adj_ += adjsigma;
        }
    };

    /**
     * This is a subclass of the vari class for precomputed
     * gradients of cov_exp_quad, specialized for a
     * double valued amplitude, sigma.
     *
     * The class stores the double values for the distance
     * matrix, pointers to the varis for the covariance 
     * matrix, along with a pointer to the vari for l.
     *
     * @tparam T_x type of std::vector of elements
     * @tparam T_l type of length scale
     */
    template <typename T_x, typename T_l>
    class cov_exp_quad_vari<T_x, double, T_l> : public vari {
      public:
        int size_;
        double l_d_;
        double sigma_d_;
        double sigma_sq_d_;
        double* dist_;
        vari* l_vari_;
        vari** cov_;

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
                          const double sigma,
                          const T_l& l)
          : vari(0.0),
          size_(x.size()),
          l_d_(value_of(l)), sigma_d_(value_of(sigma)),
          sigma_sq_d_(std::pow(sigma_d_, 2)),
          dist_(ChainableStack::memalloc_.alloc_array<double>(size_
                                                              * size_)),
          l_vari_(l.vi_),
          cov_(ChainableStack::memalloc_.alloc_array<vari*>(size_
                                                            * size_)) {
            size_t pos = 0;
            double inv_half_sq_l_d = 0.5 / (std::pow(l_d_, 2));
            for (size_t j = 0; j < static_cast<size_t>(size_); ++j)
              for (size_t i = 0; i < static_cast<size_t>(size_); ++i) {
                dist_[pos] = squared_distance(x[i], x[j])
                             * inv_half_sq_l_d;
                cov_[pos] = new vari((i == j) ? sigma_sq_d_
                                     : sigma_sq_d_ * exp(-dist_[pos]), false);
                ++pos;
              }
          }

        virtual void chain() {
          using Eigen::MatrixXd;
          using Eigen::ArrayXXd;
          using Eigen::Map;
          double adjl;
          ArrayXXd adj_cov(size_, size_);
          ArrayXXd cov(size_, size_);

          for (size_t i = 0; i < static_cast<size_t>(adj_cov.size()); ++i) {
            adj_cov(i) = cov_[i]->adj_;
            cov(i) = cov_[i]->val_;
          }
          adjl = (adj_cov * Map<ArrayXXd>(dist_, size_, size_)
            * cov).sum() * 2 / l_d_;
          l_vari_->adj_ += adjl;
        }
    };

    /**
     * This is a subclass of the vari class for precomputed
     * gradients of cov_exp_quad, specialized for a
     * double valued length scale, l.
     *
     * The class stores the double values for the distance
     * matrix, pointers to the varis for the covariance 
     * matrix, along with a pointer to the vari for sigma.
     *
     * @tparam T_x type of std::vector of elements
     * @tparam T_sigma type of sigma
     */
    template <typename T_x, typename T_sigma>
    class cov_exp_quad_vari<T_x, T_sigma, double> : public vari {
      public:
        int size_;
        double l_d_;
        double sigma_d_;
        double sigma_sq_d_;
        double* dist_;
        vari* sigma_vari_;
        vari** cov_;

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
                          const double l)
          : vari(0.0),
          size_(x.size()),
          l_d_(value_of(l)), sigma_d_(value_of(sigma)),
          sigma_sq_d_(std::pow(sigma_d_, 2)),
          dist_(ChainableStack::memalloc_.alloc_array<double>(size_
                                                                  * size_)),
          sigma_vari_(sigma.vi_),
          cov_(ChainableStack::memalloc_.alloc_array<vari*>(size_
                                                            * size_)) {
            size_t pos = 0;
            double inv_half_sq_l_d = 0.5 / (std::pow(l_d_, 2));
            for (size_t j = 0; j < static_cast<size_t>(size_); ++j)
              for (size_t i = 0; i < static_cast<size_t>(size_); ++i) {
                dist_[pos] = squared_distance(x[i], x[j])
                             * inv_half_sq_l_d;
                cov_[pos] = new vari((i == j) ? sigma_sq_d_
                                     : sigma_sq_d_ * exp(-dist_[pos]), false);
                ++pos;
              }
          }

        virtual void chain() {
          using Eigen::MatrixXd;
          using Eigen::ArrayXXd;
          using Eigen::Map;
          double adjsigma;
          ArrayXXd adj_cov(size_, size_);
          ArrayXXd cov(size_, size_);

          for (size_t i = 0; i < static_cast<size_t>(adj_cov.size()); ++i) {
            adj_cov(i) = cov_[i]->adj_;
            cov(i) = cov_[i]->val_;
          }
          adjsigma = 2 / sigma_d_ * (adj_cov * cov).sum();
          sigma_vari_->adj_ += adjsigma;
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
    template <typename T_x, typename T_sigma, typename T_l>
    inline typename
    boost::enable_if_c<(boost::is_same<T_sigma, var>::value
                       || boost::is_same<T_l, var>::value)
                       & boost::is_same<typename scalar_type<T_x>::type,
                                        double>::value,
                       Eigen::Matrix<var, -1, -1> >::type
    cov_exp_quad(const std::vector<T_x>& x,
                 const T_sigma& sigma,
                 const T_l& l) {
      check_positive("cov_exp_quad", "sigma", sigma);
      check_positive("cov_exp_quad", "l", l);
      for (size_t n = 0; n < x.size(); n++)
        check_not_nan("cov_exp_quad", "x", x[n]);

      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
        cov(x.size(), x.size());
      if (x.size() == 0)
        return cov;

      cov_exp_quad_vari<T_x, T_sigma, T_l> *baseVari
        = new cov_exp_quad_vari<T_x, T_sigma, T_l>(x, sigma, l);

      for (size_t i = 0; i < static_cast<size_t>(cov.size()); ++i)
        cov.coeffRef(i).vi_ = baseVari->cov_[i];
      return cov;
    }
  }
}
#endif
