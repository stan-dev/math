#ifndef STAN_MATH_PRIM_FUN_MAKE_NU_HPP
#define STAN_MATH_PRIM_FUN_MAKE_NU_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the degrees of freedom for the t distribution that
 * corresponds to the shape parameter in the LKJ distribution.
 *
 * @tparam T scalar type
 * @param eta LKJ distribution parameter in (0, inf)
 * @param K number of variables in correlation matrix
 */
template <typename T>
Eigen::Array<T, Eigen::Dynamic, 1> make_nu(const T& eta, size_t K) {
  using size_type = index_type_t<Eigen::Matrix<T, Eigen::Dynamic, 1>>;

  // Best (1978) implies nu = 2 * alpha for the dof in a t
  // distribution that generates a beta variate on (-1, 1)
  Eigen::Array<T, Eigen::Dynamic, 1> nu(K * (K - 1) / 2);
  T alpha = eta + 0.5 * (K - 2.0);  // from Lewandowski et. al.
  T alpha2 = 2.0 * alpha;
  for (size_type j = 0; j < (K - 1); ++j) {
    nu(j) = alpha2;
  }
  size_t counter = K - 1;
  for (size_type i = 1; i < (K - 1); ++i) {
    alpha -= 0.5;
    alpha2 = 2.0 * alpha;
    for (size_type j = i + 1; j < K; ++j, ++counter) {
      nu(counter) = alpha2;
    }
  }
  return nu;
}

}  // namespace math
}  // namespace stan

#endif
