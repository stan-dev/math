#ifndef STAN_MATH_PRIM_FUN_FACTOR_U_HPP
#define STAN_MATH_PRIM_FUN_FACTOR_U_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/atanh.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

/**
 * This function is intended to make starting values, given a unit
 * upper-triangular matrix U such that U'DU is a correlation matrix
 *
 * @tparam T type of elements in the matrix
 * @param U Sigma matrix
 * @param CPCs fill this unbounded
 */
template <typename T_U, typename T_CPCs, require_eigen_t<T_U>* = nullptr,
          require_eigen_vector_t<T_CPCs>* = nullptr,
          require_vt_same<T_U, T_CPCs>* = nullptr>
void factor_U(const T_U& U, T_CPCs&& CPCs) {
  size_t K = U.rows();
  size_t position = 0;
  size_t pull = K - 1;

  const Eigen::Ref<const plain_type_t<T_U>>& U_ref = U;

  if (K == 2) {
    CPCs(0) = atanh(U_ref(0, 1));
    return;
  }

  Eigen::Array<value_type_t<T_U>, 1, Eigen::Dynamic> temp
      = U_ref.row(0).tail(pull);

  CPCs.head(pull) = temp;

  Eigen::Array<value_type_t<T_U>, Eigen::Dynamic, 1> acc(K);
  acc(0) = -0.0;
  acc.tail(pull) = 1.0 - temp.square();
  for (size_t i = 1; i < (K - 1); i++) {
    position += pull;
    pull--;
    temp = U_ref.row(i).tail(pull);
    temp /= sqrt(acc.tail(pull) / acc(i));
    CPCs.segment(position, pull) = temp;
    acc.tail(pull) *= 1.0 - temp.square();
  }
  CPCs = 0.5 * ((1.0 + CPCs) / (1.0 - CPCs)).log();  // now unbounded
}

}  // namespace math
}  // namespace stan

#endif
