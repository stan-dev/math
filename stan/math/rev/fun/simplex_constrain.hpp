#ifndef STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/adj_jac_apply.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename T>
class simplex_constrain_op {
 public:
  adj_arg_t<T> diag_;
  adj_arg_t<T> z_;
  explicit simplex_constrain_op(const T& x)
      : diag_(setup_adj_arg<T>(x.size())), z_(setup_adj_arg<T>(x.size())) {}
  /**
   * Return the simplex corresponding to the specified free vector.
   * A simplex is a vector containing values greater than or equal
   * to 0 that sum to 1.  A vector with (K-1) unconstrained values
   * will produce a simplex of size K.
   *
   * The transform is based on a centered stick-breaking process.
   *
   * @tparam ColVec type of the vector
   * @param y Free vector input of dimensionality K - 1
   * @return Simplex of dimensionality K
   */
  template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
  auto operator()(ColVec&& y) {
    ref_type_t<ColVec> y_ref(std::forward<ColVec>(y));
    Eigen::Matrix<double, Eigen::Dynamic, 1> x(y_ref.size() + 1);
    double stick_len(1.0);
    for (int k = 0; k < y_ref.size(); ++k) {
      double log_N_minus_k = std::log(y_ref.size() - k);
      z_(k) = inv_logit(y_ref(k) - log_N_minus_k);
      diag_(k) = stick_len * z_(k) * inv_logit(log_N_minus_k - y_ref(k));
      x(k) = stick_len * z_(k);
      stick_len -= x(k);
    }
    x(y_ref.size()) = stick_len;
    return x;
  }

  /**
   * Compute the result of multiply the transpose of the adjoint vector times
   * the Jacobian of the simplex_constrain operator.
   *
   * @tparam ColVec type of the vector
   * @param adj Eigen::VectorXd of adjoints at the output of the softmax
   * @return Eigen::VectorXd of adjoints propagated through softmax operation
   */
  template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
  auto multiply_adjoint_jacobian(ColVec&& adj) const {
    ref_type_t<ColVec> adj_ref(std::forward<ColVec>(adj));
    const auto z_size = z_.size();
    Eigen::VectorXd adj_times_jac(z_size);
    double acc = adj_ref(z_size);

    if (z_size > 0) {
      adj_times_jac(z_size - 1)
          = diag_(z_size - 1) * (adj_ref(z_size - 1) - acc);
      for (int n = z_size - 1; --n >= 0;) {
        acc = adj_ref(n + 1) * z_(n + 1) + (1 - z_(n + 1)) * acc;
        adj_times_jac(n) = diag_(n) * (adj_ref(n) - acc);
      }
    }
    return std::make_tuple(adj_times_jac);
  }
};
}  // namespace internal

/**
 * Return the simplex corresponding to the specified free vector.
 * A simplex is a vector containing values greater than or equal
 * to 0 that sum to 1.  A vector with (K-1) unconstrained values
 * will produce a simplex of size K.
 *
 * The transform is based on a centered stick-breaking process.
 * @tparam ColVec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @return Simplex of dimensionality K
 */
template <typename ColVec,
          require_eigen_col_vector_vt<is_var, ColVec>* = nullptr>
inline auto simplex_constrain(ColVec&& y) {
  return adj_jac_apply<internal::simplex_constrain_op<ColVec>>(
      std::forward<ColVec>(y));
}

}  // namespace math
}  // namespace stan
#endif
