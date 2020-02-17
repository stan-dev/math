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
template <typename T, typename = require_eigen_t<T>>
inline plain_type_t<T> matrix_exp(const T& A_in) {
  using std::exp;
  const auto& A = A_in.eval();
  check_square("matrix_exp", "input matrix", A);
  if (T::RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1) {
    plain_type_t<T> res;
    res << exp(A(0));
    return res;
  }
  if (A_in.size() == 0) {
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
