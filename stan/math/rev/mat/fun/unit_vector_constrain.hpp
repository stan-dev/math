#ifndef STAN_MATH_PRIM_MAT_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_MAT_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/dot_self.hpp>
#include <stan/math/prim/mat/err/check_vector.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/dot_self.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
class unit_vector_elt_vari : public vari {
 private:
  vari** y_;
  const double* unit_vector_y_;
  const int size_;
  const int idx_;
  const double norm_;

 public:
  unit_vector_elt_vari(double val, vari** y, const double* unit_vector_y,
                       int size, int idx, double norm)
      : vari(val),
        y_(y),
        unit_vector_y_(unit_vector_y),
        size_(size),
        idx_(idx),
        norm_(norm) {}
  void chain() {
    const double cubed_norm = std::pow(norm_, 3);
    for (int m = 0; m < size_; ++m) {
      y_[m]->adj_
          -= adj_ * unit_vector_y_[m] * unit_vector_y_[idx_] / cubed_norm;
      if (m == idx_)
        y_[m]->adj_ += adj_ / norm_;
    }
  }
};
}  // namespace internal

/**
 * Return the unit length vector corresponding to the free vector y.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 * @tparam T Scalar type.
 **/
template <int R, int C>
Eigen::Matrix<var, R, C> unit_vector_constrain(
    const Eigen::Matrix<var, R, C>& y) {
  check_vector("unit_vector", "y", y);
  check_nonzero_size("unit_vector", "y", y);
  vector_d y_d = y.val();

  vari** y_vi_array = reinterpret_cast<vari**>(
      ChainableStack::instance_->memalloc_.alloc(sizeof(vari*) * y.size()));
  double* unit_vector_y_d_array = reinterpret_cast<double*>(
      ChainableStack::instance_->memalloc_.alloc(sizeof(double) * y_d.size()));

  Eigen::Map<vector_vi>(y_vi_array, y.size()) = y.vi();
  const double norm = y_d.norm();

  check_positive_finite("unit_vector", "norm", norm);
  Eigen::Map<vector_d> unit_vecd(unit_vector_y_d_array, y.size());
  unit_vecd = y_d / norm;

  Eigen::Matrix<var, R, C> unit_vector_y(y.size());
  for (int k = 0; k < y.size(); ++k)
    unit_vector_y.coeffRef(k) = var(new internal::unit_vector_elt_vari(
        unit_vecd[k], y_vi_array, unit_vector_y_d_array, y.size(), k, norm));
  return unit_vector_y;
}

/**
 * Return the unit length vector corresponding to the free vector y.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 * @param lp Log probability reference to increment.
 * @tparam T Scalar type.
 **/
template <int R, int C>
Eigen::Matrix<var, R, C> unit_vector_constrain(
    const Eigen::Matrix<var, R, C>& y, var& lp) {
  Eigen::Matrix<var, R, C> x = unit_vector_constrain(y);
  lp -= 0.5 * dot_self(y);
  return x;
}

}  // namespace math
}  // namespace stan
#endif
