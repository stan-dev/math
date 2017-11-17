#ifndef STAN_MATH_REV_MAT_FUN_COV_PERIODIC_HPP
#define STAN_MATH_REV_MAT_FUN_COV_PERIODIC_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <vector>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * This is a subclass of the vari class for precomputed
     * gradients of cov_periodic.
     *
     * The class stores the double values for the distance
     * matrix, pointers to the varis for the covariance
     * matrix, along with a pointer to the vari for sigma,
     * the vari for l and the vari for p.
     *
     * @tparam T_x type of std::vector of elements
     * @tparam T_sigma type of sigma
     * @tparam T_l type of length scale
     * @tparam T_p type of period
     */
    template <typename T_x, typename T_sigma, typename T_l, typename T_p>
    class cov_periodic_vari : public vari {
    public:
      const size_t size_;
      const size_t size_ltri_;
      const double l_d_;
      const double sigma_d_;
      const double p_d_;
      const double sigma_sq_d_;
      double* dist_;
      double* sin_dist_;
      double* cos_dist_;
      double* sin_dist_sq_;
      vari* l_vari_;
      vari* sigma_vari_;
      vari* p_vari_;
      vari** cov_lower_;
      vari** cov_diag_;

      /**
       * Constructor for cov_periodic.
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
       * @param p period
       */
      cov_periodic_vari(const std::vector<T_x>& x,
                        const T_sigma& sigma,
                        const T_l& l,
                        const T_p& p
						)
        : vari(0.0),
          size_(x.size()),
          size_ltri_(size_ * (size_ - 1) / 2),
          l_d_(value_of(l)), sigma_d_(value_of(sigma)),
		  p_d_(value_of(p)),
          sigma_sq_d_(sigma_d_ * sigma_d_),
		  dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
		  sin_dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
		  cos_dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
		  sin_dist_sq_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
          l_vari_(l.vi_), sigma_vari_(sigma.vi_),
		  p_vari_(p.vi_),
          cov_lower_(ChainableStack::memalloc_.alloc_array<vari*>(size_ltri_)),
          cov_diag_(ChainableStack::memalloc_.alloc_array<vari*>(size_)) {
        double neg_two_inv_l_sq = -2.0 / (l_d_ * l_d_);
        double pi_div_p = M_PI / p_d_;

        size_t pos = 0;
        for (size_t j = 0; j < size_ - 1; ++j) {
          for (size_t i = j + 1; i < size_; ++i) {
        	double dist = distance(x[i], x[j]);
        	double sin_dist = sin(pi_div_p * dist);
            double sin_dist_sq = square(sin_dist);
            dist_[pos] = dist;
            sin_dist_[pos] = sin_dist;
        	cos_dist_[pos] = cos(pi_div_p * dist);
            sin_dist_sq_[pos] = sin_dist_sq;
            cov_lower_[pos] = new vari(sigma_sq_d_ * std::exp(sin_dist_sq
                                       * neg_two_inv_l_sq), false);
            ++pos;
          }
        }
        for (size_t i = 0; i < size_; ++i)
          cov_diag_[i] = new vari(sigma_sq_d_, false);
      }

      virtual void chain() {
        double adjl = 0;
        double adjsigma = 0;
        double adjp = 0;

        for (size_t i = 0; i < size_ltri_; ++i) {
          vari* el_low = cov_lower_[i];
          double prod_add = el_low->adj_ * el_low->val_;
          adjl += prod_add * sin_dist_sq_[i];
          adjsigma += prod_add;
          adjp += prod_add * sin_dist_[i] * cos_dist_[i] * dist_[i];
        }
        for (size_t i = 0; i < size_; ++i) {
          vari* el = cov_diag_[i];
          adjsigma += el->adj_ * el->val_;
        }
        double l_d_sq = l_d_ * l_d_;
        l_vari_->adj_ +=  adjl * 4 / (l_d_sq * l_d_);
        sigma_vari_->adj_ += adjsigma * 2 / sigma_d_;
        p_vari_->adj_ +=  adjp * 4 * M_PI / l_d_sq / (p_d_ * p_d_);
      }
    };

    /**
     * This is a subclass of the vari class for precomputed
     * gradients of cov_periodic.
     *
     * The class stores the double values for the distance
     * matrix, pointers to the varis for the covariance
     * matrix, along with a pointer to the vari for sigma,
     * the vari for l and the vari for p.
     *
     * @tparam T_x type of std::vector of elements
     * @tparam T_l type of length scale
     * @tparam T_p type of period
     */
    template <typename T_x, typename T_l, typename T_p>
    class cov_periodic_vari<T_x, double, T_l, T_p> : public vari {
    public:
      const size_t size_;
      const size_t size_ltri_;
      const double l_d_;
      const double sigma_d_;
      const double p_d_;
      const double sigma_sq_d_;
      double* dist_;
      double* sin_dist_;
      double* cos_dist_;
      double* sin_dist_sq_;
      vari* l_vari_;
      vari* p_vari_;
      vari** cov_lower_;
      vari** cov_diag_;

      /**
       * Constructor for cov_periodic.
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
       * @param p period
       */
      cov_periodic_vari(const std::vector<T_x>& x,
                        double sigma,
                        const T_l& l,
						const T_p& p)
        : vari(0.0),
          size_(x.size()),
          size_ltri_(size_ * (size_ - 1) / 2),
          l_d_(value_of(l)), sigma_d_(value_of(sigma)),
		  p_d_(value_of(p)),
          sigma_sq_d_(sigma_d_ * sigma_d_),
          dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
		  sin_dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
		  cos_dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
		  sin_dist_sq_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
          l_vari_(l.vi_),
		  p_vari_(p.vi_),
          cov_lower_(ChainableStack::memalloc_.alloc_array<vari*>(size_ltri_)),
          cov_diag_(ChainableStack::memalloc_.alloc_array<vari*>(size_)) {

        double neg_two_inv_l_sq = -2.0 / (l_d_ * l_d_);
        double pi_div_p = M_PI / p_d_;

        size_t pos = 0;
        for (size_t j = 0; j < size_ - 1; ++j) {
          for (size_t i = j + 1; i < size_; ++i) {
          	double dist = distance(x[i], x[j]);
          	double sin_dist = sin(pi_div_p * dist);
            double sin_dist_sq = square(sin_dist);
            dist_[pos] = dist;
            sin_dist_[pos] = sin_dist;
          	cos_dist_[pos] = cos(pi_div_p * dist);
            sin_dist_sq_[pos] = sin_dist_sq;
            cov_lower_[pos] = new vari(sigma_sq_d_ * std::exp(sin_dist_sq
                                         * neg_two_inv_l_sq), false);
            ++pos;
          }
        }
        for (size_t i = 0; i < size_; ++i)
          cov_diag_[i] = new vari(sigma_sq_d_, false);
      }

      virtual void chain() {
        double adjl = 0;
        double adjp = 0;

        for (size_t i = 0; i < size_ltri_; ++i) {
          vari* el_low = cov_lower_[i];
          double prod_add = el_low->adj_ * el_low->val_;
          adjl += prod_add * sin_dist_sq_[i];
          adjp += prod_add * sin_dist_[i] * cos_dist_[i] * dist_[i];
        }
        double l_d_sq = l_d_ * l_d_;
        l_vari_->adj_ +=  adjl * 4  / (l_d_sq * l_d_);
        p_vari_->adj_ +=  adjp * 4 * M_PI / l_d_sq / (p_d_ * p_d_);
      }
    };

    /**
     * Returns a periodic kernel.
     *
     * @param x std::vector input that can be used in square distance
     *    Assumes each element of x is the same size
     * @param sigma standard deviation
     * @param l length scale
     * @param p period
     * @return squared distance
     * @throw std::domain_error if sigma <= 0, l <= 0, p <= 0, or
     *   x is nan or infinite
     */
    template <typename T_x>
    inline typename
    boost::enable_if_c<boost::is_same<typename scalar_type<T_x>::type,
                                      double>::value,
                       Eigen::Matrix<var, -1, -1> >::type
      cov_periodic(const std::vector<T_x>& x,
                   const var& sigma,
                   const var& l,
				   const var& p) {
      check_positive("cov_periodic", "marginal variance", sigma);
      check_positive("cov_periodic", "length-scale", l);
      check_positive("cov_periodic", "period", p);
      size_t x_size = x.size();
      for (size_t i = 0; i < x_size; ++i)
        check_not_nan("cov_periodic", "x", x[i]);

      Eigen::Matrix<var, -1, -1>
        cov(x_size, x_size);
      if (x_size == 0)
        return cov;

      cov_periodic_vari<T_x, var, var, var> *baseVari
        = new cov_periodic_vari<T_x, var, var, var>(x, sigma, l, p);

      size_t pos = 0;
      for (size_t j = 0; j < x_size - 1; ++j) {
        for (size_t i = (j + 1); i < x_size; ++i) {
          cov.coeffRef(i, j).vi_ = baseVari->cov_lower_[pos];
          cov.coeffRef(j, i).vi_ = cov.coeffRef(i, j).vi_;
          ++pos;
        }
        cov.coeffRef(j, j).vi_ = baseVari->cov_diag_[j];
      }
      cov.coeffRef(x_size - 1, x_size - 1).vi_
        = baseVari->cov_diag_[x_size - 1];
      return cov;
    }

    /**
     * Returns a periodic kernel.
     *
     * @param x std::vector input that can be used in square distance
     *    Assumes each element of x is the same size
     * @param sigma standard deviation
     * @param l length scale
     * @param p period
     * @return squared distance
     * @throw std::domain_error if sigma <= 0, l <= 0, p <= 0, or
     *   x is nan or infinite
     */
    template <typename T_x>
    inline typename
    boost::enable_if_c<boost::is_same<typename scalar_type<T_x>::type,
                                      double>::value,
                       Eigen::Matrix<var, -1, -1> >::type
      cov_periodic(const std::vector<T_x>& x,
                   double sigma,
                   const var& l,
				   const var& p) {

      check_positive("cov_periodic", "marginal variance", sigma);
      check_positive("cov_periodic", "length-scale", l);
      check_positive("cov_periodic", "period", p);
      size_t x_size = x.size();
      for (size_t i = 0; i < x_size; ++i)
        check_not_nan("cov_periodic", "x", x[i]);

      Eigen::Matrix<var, -1, -1>
        cov(x_size, x_size);
      if (x_size == 0)
        return cov;

      cov_periodic_vari<T_x, double, var, var> *baseVari
        = new cov_periodic_vari<T_x, double, var, var>(x, sigma, l, p);

      size_t pos = 0;
      for (size_t j = 0; j < x_size - 1; ++j) {
        for (size_t i = (j + 1); i < x_size; ++i) {
          cov.coeffRef(i, j).vi_ = baseVari->cov_lower_[pos];
          cov.coeffRef(j, i).vi_ = cov.coeffRef(i, j).vi_;
          ++pos;
        }
        cov.coeffRef(j, j).vi_ = baseVari->cov_diag_[j];
      }
      cov.coeffRef(x_size - 1, x_size - 1).vi_
        = baseVari->cov_diag_[x_size - 1];
      return cov;
    }

  }
}
#endif
