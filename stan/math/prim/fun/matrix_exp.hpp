#ifndef STAN_MATH_PRIM_FUN_MATRIX_EXP_HPP
#define STAN_MATH_PRIM_FUN_MATRIX_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/matrix_exp_pade.hpp>
#include <stan/math/prim/fun/matrix_exp_2x2.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the matrix exponential of the input
 * matrix.
 *
 * @tparam T type of the matrix
 * @param[in] A_in Matrix to exponentiate.
 * @return Matrix exponential, dynamically-sized.
 * @throw <code>std::invalid_argument</code> if the input matrix
 * is not square.
 */
template <typename EigenMat, typename = require_eigen_t<EigenMat>>
inline plain_type_t<EigenMat> matrix_exp(EigenMat&& A_in) {
  decltype(auto) A = to_ref(std::forward<EigenMat>(A_in));
  using T = std::decay_t<EigenMat>;
  check_square("matrix_exp", "input matrix", A);
  if constexpr (T::RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1) {
    plain_type_t<T> res;
    res << exp(A(0));
    return res;
  }
  if (A.size() == 0) {
    return {};
  }
  return (A.cols() == 2
          && square(value_of(A(0, 0)) - value_of(A(1, 1)))
                     + 4 * value_of(A(0, 1)) * value_of(A(1, 0))
                 > 0)
             ? matrix_exp_2x2(A)
             : matrix_exp_pade(A);
}

}  // namespace math
}  // namespace stan

#endif
