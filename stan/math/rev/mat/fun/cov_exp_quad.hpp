#ifndef STAN_MATH_REV_MAT_FUN_COV_EXP_QUAD_HPP
#define STAN_MATH_REV_MAT_FUN_COV_EXP_QUAD_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>
#include <stan/math/prim/scal/fun/exp.hpp>
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
 * gradients of cov_exp_quad.
 *
 * The class stores the double values for the distance
 * matrix, pointers to the varis for the covariance
 * matrix, along with a pointer to the vari for sigma,
 * and the vari for l.
 *
 * @deprecated use <code>gp_exp_quad_cov_vari</code>
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_sigma type of sigma
 * @tparam T_l type of length scale
 */
template <typename T_x, typename T_sigma, typename T_l>
class cov_exp_quad_vari : public vari {
 public:
  const size_t size_;
  const size_t size_ltri_;
  const double l_d_;
  const double sigma_d_;
  const double sigma_sq_d_;
  double* dist_;
  vari* l_vari_;
  vari* sigma_vari_;
  vari** cov_lower_;
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
   * @deprecated use <code>gp_exp_quad_cov_vari</code>
   *
   * @param x std::vector input that can be used in square distance
   *    Assumes each element of x is the same size
   * @param sigma standard deviation
   * @param l length scale
   */
  cov_exp_quad_vari(const std::vector<T_x>& x, const T_sigma& sigma,
                    const T_l& l)
      : vari(0.0),
        size_(x.size()),
        size_ltri_(size_ * (size_ - 1) / 2),
        l_d_(value_of(l)),
        sigma_d_(value_of(sigma)),
        sigma_sq_d_(sigma_d_ * sigma_d_),
        dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
        l_vari_(l.vi_),
        sigma_vari_(sigma.vi_),
        cov_lower_(ChainableStack::memalloc_.alloc_array<vari*>(size_ltri_)),
        cov_diag_(ChainableStack::memalloc_.alloc_array<vari*>(size_)) {
    double inv_half_sq_l_d = 0.5 / (l_d_ * l_d_);
    size_t pos = 0;
    for (size_t j = 0; j < size_ - 1; ++j) {
      for (size_t i = j + 1; i < size_; ++i) {
        double dist_sq = squared_distance(x[i], x[j]);
        dist_[pos] = dist_sq;
        cov_lower_[pos] = new vari(
            sigma_sq_d_ * std::exp(-dist_sq * inv_half_sq_l_d), false);
        ++pos;
      }
    }
    for (size_t i = 0; i < size_; ++i)
      cov_diag_[i] = new vari(sigma_sq_d_, false);
  }

  virtual void chain() {
    double adjl = 0;
    double adjsigma = 0;

    for (size_t i = 0; i < size_ltri_; ++i) {
      vari* el_low = cov_lower_[i];
      double prod_add = el_low->adj_ * el_low->val_;
      adjl += prod_add * dist_[i];
      adjsigma += prod_add;
    }
    for (size_t i = 0; i < size_; ++i) {
      vari* el = cov_diag_[i];
      adjsigma += el->adj_ * el->val_;
    }
    l_vari_->adj_ += adjl / (l_d_ * l_d_ * l_d_);
    sigma_vari_->adj_ += adjsigma * 2 / sigma_d_;
  }
};

/**
 * This is a subclass of the vari class for precomputed
 * gradients of cov_exp_quad.
 *
 * The class stores the double values for the distance
 * matrix, pointers to the varis for the covariance
 * matrix, along with a pointer to the vari for sigma,
 * and the vari for l.
 *
 * @deprecated use <code>gp_exp_quad_cov_vari</code>
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_l type of length scale
 */
template <typename T_x, typename T_l>
class cov_exp_quad_vari<T_x, double, T_l> : public vari {
 public:
  const size_t size_;
  const size_t size_ltri_;
  const double l_d_;
  const double sigma_d_;
  const double sigma_sq_d_;
  double* dist_;
  vari* l_vari_;
  vari** cov_lower_;
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
   * @deprecated use <code>gp_exp_quad_cov_vari</code>
   *
   * @param x std::vector input that can be used in square distance
   *    Assumes each element of x is the same size
   * @param sigma standard deviation
   * @param l length scale
   */
  cov_exp_quad_vari(const std::vector<T_x>& x, double sigma, const T_l& l)
      : vari(0.0),
        size_(x.size()),
        size_ltri_(size_ * (size_ - 1) / 2),
        l_d_(value_of(l)),
        sigma_d_(value_of(sigma)),
        sigma_sq_d_(sigma_d_ * sigma_d_),
        dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
        l_vari_(l.vi_),
        cov_lower_(ChainableStack::memalloc_.alloc_array<vari*>(size_ltri_)),
        cov_diag_(ChainableStack::memalloc_.alloc_array<vari*>(size_)) {
    double inv_half_sq_l_d = 0.5 / (l_d_ * l_d_);
    size_t pos = 0;
    for (size_t j = 0; j < size_ - 1; ++j) {
      for (size_t i = j + 1; i < size_; ++i) {
        double dist_sq = squared_distance(x[i], x[j]);
        dist_[pos] = dist_sq;
        cov_lower_[pos] = new vari(
            sigma_sq_d_ * std::exp(-dist_sq * inv_half_sq_l_d), false);
        ++pos;
      }
    }
    for (size_t i = 0; i < size_; ++i)
      cov_diag_[i] = new vari(sigma_sq_d_, false);
  }

  virtual void chain() {
    double adjl = 0;

    for (size_t i = 0; i < size_ltri_; ++i) {
      vari* el_low = cov_lower_[i];
      adjl += el_low->adj_ * el_low->val_ * dist_[i];
    }
    l_vari_->adj_ += adjl / (l_d_ * l_d_ * l_d_);
  }
};

/**
 * Returns a squared exponential kernel.
 *
 * @deprecated use <code>gp_exp_quad_cov_vari</code>
 *
 * @param x std::vector input that can be used in square distance
 *    Assumes each element of x is the same size
 * @param sigma standard deviation
 * @param l length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x>
inline typename boost::enable_if_c<
    boost::is_same<typename scalar_type<T_x>::type, double>::value,
    Eigen::Matrix<var, -1, -1> >::type
cov_exp_quad(const std::vector<T_x>& x, const var& sigma, const var& l) {
  return gp_exp_quad_cov(x, sigma, l);
}

/**
 * Returns a squared exponential kernel.
 *
 * @deprecated use <code>gp_exp_quad_cov_vari</code>
 *
 * @param x std::vector input that can be used in square distance
 *    Assumes each element of x is the same size
 * @param sigma standard deviation
 * @param l length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x>
inline typename boost::enable_if_c<
    boost::is_same<typename scalar_type<T_x>::type, double>::value,
    Eigen::Matrix<var, -1, -1> >::type
cov_exp_quad(const std::vector<T_x>& x, double sigma, const var& l) {
  return gp_exp_quad_cov(x, sigma, l);
}

}  // namespace math
}  // namespace stan
#endif
