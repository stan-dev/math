#ifndef STAN_MATH_REV_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_REV_FUN_LOG_SOFTMAX_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log_softmax.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

namespace internal {

class log_softmax_elt_vari : public vari {
 private:
  vari** alpha_;
  const double* softmax_alpha_;
  const int size_;  // array sizes
  const int idx_;   // in in softmax output

 public:
  log_softmax_elt_vari(double val, vari** alpha, const double* softmax_alpha,
                       int size, int idx)
      : vari(val),
        alpha_(alpha),
        softmax_alpha_(softmax_alpha),
        size_(size),
        idx_(idx) {}
  void chain() {
    for (int m = 0; m < size_; ++m) {
      if (m == idx_) {
        alpha_[m]->adj_ += adj_ * (1 - softmax_alpha_[m]);
      } else {
        alpha_[m]->adj_ -= adj_ * softmax_alpha_[m];
      }
    }
  }
};
}  // namespace internal

/**
 * Return the log softmax of the specified vector
 *
 * @tparam T type of input
 * @param x input
 * @return softmax of the input
 * @throw std::domain_error if the input size is 0
 */
template <typename T, require_eigen_st<is_var, T>* = nullptr>
auto log_softmax(const T& x) {
  const int a_size = x.size();

  check_nonzero_size("log_softmax", "x", x);

  const auto& x_ref = to_ref(x);

  vari** x_vi_array
      = ChainableStack::instance_->memalloc_.alloc_array<vari*>(a_size);
  Eigen::Map<vector_vi>(x_vi_array, a_size) = x_ref.vi();

  vector_d x_d = x_ref.val();

  // fold logic of math::softmax() and math::log_softmax()
  // to save computations

  vector_d diff = (x_d.array() - x_d.maxCoeff());
  vector_d softmax_x_d = diff.array().exp();
  double sum = softmax_x_d.sum();
  vector_d log_softmax_x_d = diff.array() - std::log(sum);

  // end fold
  double* softmax_x_d_array
      = ChainableStack::instance_->memalloc_.alloc_array<double>(a_size);
  Eigen::Map<vector_d>(softmax_x_d_array, a_size) = softmax_x_d.array() / sum;

  plain_type_t<T> log_softmax_x(a_size);
  for (int k = 0; k < a_size; ++k) {
    log_softmax_x(k) = var(new internal::log_softmax_elt_vari(
        log_softmax_x_d[k], x_vi_array, softmax_x_d_array, a_size, k));
  }
  return log_softmax_x;
}

/**
 * Return the log softmax of the specified vector
 *
 * @tparam T type of input
 * @param x input
 * @return softmax of the input
 * @throw std::domain_error if the input size is 0
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto log_softmax(const T& x) {
  check_nonzero_size("log_softmax", "x", x);

  const auto& theta = (x.val().array() - x.val().maxCoeff()).eval();

  return make_callback_var(
      (theta.array() - log(theta.exp().sum())).matrix(),
      [x](const auto& res) mutable {
        x.adj().noalias()
            += res.adj() - (res.adj().sum() * res.val().array().exp()).matrix();
      });
}

/**
 * Return the log softmax of the specified `std::vector` or
 * `std::vector` of containers.
 *
 * @tparam T type of input
 * @param x input
 * @return softmax of the input
 * @throw std::domain_error if the input size is 0
 */
template <typename T, require_std_vector_st<is_var, T>* = nullptr>
inline auto log_softmax(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [](const auto& alpha) { return log_softmax(alpha); });
}

}  // namespace math
}  // namespace stan
#endif
