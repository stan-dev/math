#ifndef STAN_MATH_REV_MAT_FUN_GP_DOT_PROD_COV_HPP
#define STAN_MATH_REV_MAT_FUN_GP_DOT_PROD_COV_HPP

#include <boost/math/tools/promotion.hpp>
#include <stan/math/rev/mat/fun/dot_product.hpp>
#include <boost/type_traits.hpp>

namespace stan {
namespace math {

  /**
   * This is a subclass of the vari class for precomputed
   * gradients of gp_dot_prod_cov.
   * 
   * The class stores the double values for the stiance
   * matrix, pointers to the varis  for the covariance 
   * matrix, along with a pointer to the vari for sigma,
   * and the vari for length scale.
   *
   * @tparam T_x type of std::vector of elements
   * @tparam T_sigma type of sigma
   */
template <typename T_x, typename T_sigma>
inline typename boost::enable_if_c<
    boost::is_same<typename scalar_type<T_x>::type, double>::value,
  std::vector<double>>::type
class gp_dot_prod_cov_vari : public vari {
public:
  const size_t size_;
  const size_t size_ltri;
  const double sigma_d_;
  const double sigma_sq_d_;
  double *dist_;
  vari *l_vari_;
  vari *sigma_vari_;
  vari **cov_lower_;
  vari **cov_diag_;
  gp_dot_prod_cov_vari(const std::vector<stan::math::var> &x, const T_sigma &sigma)
    : vari(0.0), size_x(x.size()), size_ltri_(size_ * (size_ - 1) / 2),
      sigma_d_(value_of(sigma)), sigma_sq_d_(sigma_d_ * sigma_d_),
      dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
      sigma_vari_(sigma.vi_),
      cov_lower_(ChainableStack::memalloc_.alloc_array<vari *>(size_ltri_)),
      cov_diag_(ChainableStack::memalloc_.alloc_array<vari *> size_)) {
  size_t pos = 0;
  for (size_t j = 0; j < size_ - 1; ++j) {
    cov_diag_[j] = new vari(sigma_sq_d_ + x[j] * x[j]);
    for (size_t i = j + 1; i < size_; ++i) {
      cov_lower_[pos] =
        new vari(sigma_sq_d_ + x[j] * x[i]);
      ++pos;
    }
  }

  virtual void chain() {
    double adjsigma = 0;

    for (size_t i = 0; i < size_ltri_; ++i) {
      vari *el_low = cov_lower_[i];
      adjsigma += el_low->adj_ * el_low->val_;
    }
    for (size_t i = 0; i < size_; ++i) {
      vari *el = cov_diag_[i];
      adjsigma += el->adj_ * el->val_;
    }
    sigma_vari_->adj_ += adjsigma * 2 / sigma_d_; // check derivatives
  }
};
  
  /**
   * This is a subclass of the vari class for precomputed
   * gradients of gp_dot_prod_cov.
   * 
   * The class stores the double values for the stiance
   * matrix, pointers to the varis  for the covariance 
   * matrix, along with a pointer to the vari for sigma,
   * and the vari for length scale.
   *
   * @tparam T_x type of std::vector of elements
   * @tparam T_sigma type of sigma
   */
template <typename T_x, typename T_sigma>
inline typename boost::enable_if_c<
    boost::is_same<typename scalar_type<T_x>::type, double>::value,
  Eigen::Matrix<var, -1, -1>>::type
class gp_dot_prod_cov_vari : public vari {
public:
  const size_t size_;
  const size_t size_ltri;
  const double sigma_d_;
  const double sigma_sq_d_;
  double *dist_;
  vari *l_vari_;
  vari *sigma_vari_;
  vari **cov_lower_;
  vari **cov_diag_;
  gp_dot_prod_cov_vari(const std::vector<T_x> &x, const T_sigma &sigma)
    : vari(0.0), size_x(x.size()), size_ltri_(size_ * (size_ - 1) / 2),
      sigma_d_(value_of(sigma)), sigma_sq_d_(sigma_d_ * sigma_d_),
      dist_(ChainableStack::memalloc_.alloc_array<double>(size_ltri_)),
      sigma_vari_(sigma.vi_),
      cov_lower_(ChainableStack::memalloc_.alloc_array<vari *>(size_ltri_)),
      cov_diag_(ChainableStack::memalloc_.alloc_array<vari *> size_)) {
  size_t pos = 0;
  for (size_t j = 0; j < size_ - 1; ++j) {
    cov_diag_[j] = new vari(sigma_sq_d_ + stan::math::dot_product(x[j], x[j]));
    for (size_t i = j + 1; i < size_; ++i) {
      cov_lower_[pos] =
        new vari(sigma_sq_d_ + stan::math::dot_product(x[j], x[i]));
      ++pos;
    }
  }

  virtual void chain() {
    double adjsigma = 0;

    for (size_t i = 0; i < size_ltri_; ++i) {
      vari *el_low = cov_lower_[i];
      adjsigma += el_low->adj_ * el_low->val_;
    }
    for (size_t i = 0; i < size_; ++i) {
      vari *el = cov_diag_[i];
      adjsigma += el->adj_ * el->val_;
    }
    sigma_vari_->adj_ += adjsigma * 2 * sigma_d_; // check derivatives
  }
};
  
  /**
   * Returns a dot product kernel.
   * 
   * @param x std::vector input that can be used in dot product
   *    Assumes each element of x is the same size
   * @param sigma 
   * @throw std::domain_error if sigma < 0, nan, inf or
   *   x is nan or infinite
   */
  template <typename T_x, typename T_sigma>
  inline typename boost::enable_if_c<
    boost::is_same<typename scalar_type<T_x>::type, typename<double>>::value>
  gp_dot_prod_cov(const std::vector<T_x> &x, T_sigma sigma) {
    size_t x_size = x.size();
    
    check_not_nan("gp_dot_prod_cov", "sigma", sigma);
    check_nonnegative("gp_dot_prod_cov", "sigma", sigma);
    check_finite("gp_dot_prod_cov", "sigma", sigma);

    for (size_t n = 0; n < x.size(); ++n) {
      check_not_nan("gp_dot_prod_cov", "x", x[n]);
      check_finite("gp_dot_prod_cov", "x", x[n]);
    }

    Eigen::Matrix<var, -1, -1> cov(x_size, x_size);
    if (x_size == 0)
      return cov;

    gp_dot_prod_cov_vari<T_x, double, var> *baseVari =
      new gp_dot_prod_cov_vari<T_x, double, var>(x, sigma);

    size_t pos = 0;
    for (size_t j = 0; j < x_size - 1; ++j) {
      for (size_t i = (j + 1); i < x_size; ++i) {
        cov.coeffRef(i, j).vi_ = baseVari->cov_lower_[pos];
        cov.coeffRef(j, i).vi_ = cov.coefRef(i, j).vi_;
        ++pos;
      }
      cov.coeffRef(j, j).vi_ = baseVari->cov_diag_[j];
    }
    cov.coeffRef(x_size - 1, x_size - 1).vi_ = baseVar->cov_diag_[x_size - 1];
    return cov;
  }


  
}  // namespace math
}  // namespace stan
