#ifndef STAN_MATH_REV_MAT_FUN_GP_EXPONENTIAL_COV_HPP
#define STAN_MATH_REV_MAT_FUN_GP_EXPONENTIAL_COV_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>
#include <stan/math/prim/scal/fun/exp.hpp>

// testing
#include <boost/utility/enable_if.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/type_traits.hpp>

namespace stan {
namespace math {

template <typename T_x, typename T_s, typename T_l>
class gp_exponential_cov_vari : public vari {
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

  gp_exponential_cov_vari(const std::vector<T_x>& x, const T_s& sigma,
                     const T_l& length_scale)
    : vari(0.0),
      size_(x.size()),
      size_ltri_(size_ * (size_ - 1) / 2),
      l_d_(value_of(length_scale)),
      sigma_d_(value_of(sigma)),
      sigma_sq_d_(sigma_d_ * sigma_d_),
      dist_(ChainableStack::context().memalloc_.alloc_array<double>(
        size_ltri_)),
      l_vari_(length_scale.vi_),
      sigma_vari_(sigma.vi_),
      cov_lower_(
            ChainableStack::context().memalloc_.alloc_array<vari*>(size_ltri_)),
      cov_diag_(
            ChainableStack::context().memalloc_.alloc_array<vari*>(size_)) {
    double neg_inv_l = -1.0 / l_d_;
    size_t pos = 0;
    for (size_t j = 0; j < size_ - 1; ++j) {
      for (size_t i = j + 1; i < size_; ++i) {
        double dist_sq = squared_distance(x[i], x[j]).val();
        dist_[pos] = dist_sq;
        cov_lower_[pos] = new vari(sigma_sq_d_ * exp(dist_sq * neg_inv_l), false);
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
      l_vari_->adj_ += adjl / (-1.0 * l_d_ * l_d_);
      sigma_vari_->adj_ += adjsigma * 2 / sigma_d_;
    }
};

template <typename T_x>
//inline typename boost::enable_if_c<
  // boost::is_same<typename scalar_type<T_x>::type, stan::math::var>::value,
  //   Eigen::Matrix<var, -1, -1> >::type
  inline typename Eigen::Matrix<typename stan::return_type<T_x>::type,
                                Eigen::Dynamic, Eigen::Dynamic>
 gp_exponential_cov(const std::vector<T_x>& x, const var& sigma, const var& l) {
  check_positive("gp_exponential_cov", "sigma", sigma);
  check_positive("gp_exponential_cov", "l", l);
  size_t x_size = x.size();
  for (size_t i = 0; i < x_size; ++i)
    check_not_nan("gp_exponential_cov", "x", x[i]);

  Eigen::Matrix<var, -1, -1> cov(x_size, x_size);
  if (x_size == 0)
    return cov;

  gp_exponential_cov_vari<T_x, var, var>* baseVari
      = new gp_exponential_cov_vari<T_x, var, var>(x, sigma, l);

  size_t pos = 0;
  for (size_t j = 0; j < x_size - 1; ++j) {
    for (size_t i = (j + 1); i < x_size; ++i) {
      cov.coeffRef(i, j).vi_ = baseVari->cov_lower_[pos];
      cov.coeffRef(j, i).vi_ = cov.coeffRef(i, j).vi_;
      ++pos;
    }
    cov.coeffRef(j, j).vi_ = baseVari->cov_diag_[j];
  }
  cov.coeffRef(x_size - 1, x_size - 1).vi_ = baseVari->cov_diag_[x_size - 1];
  return cov;
}

}  // namespace math
}  // namespace stan


#endif
