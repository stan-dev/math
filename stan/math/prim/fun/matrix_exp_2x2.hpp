#ifndef STAN_MATH_PRIM_FUN_MATRIX_EXP_2X2_HPP
#define STAN_MATH_PRIM_FUN_MATRIX_EXP_2X2_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/cosh.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sinh.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the matrix exponential of a 2x2 matrix. Reference for
 * algorithm: http://mathworld.wolfram.com/MatrixExponential.html
 * Note: algorithm only works if delta > 0;
 *
 * @tparam EigMat type of the matrix
 * @param[in] A 2x2 matrix to exponentiate.
 * @return Matrix exponential of A.
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>
matrix_exp_2x2(const EigMat& A) {
  using std::cosh;
  using std::exp;
  using std::sinh;
  using std::sqrt;

  using T = value_type_t<EigMat>;
  T a = A(0, 0), b = A(0, 1), c = A(1, 0), d = A(1, 1), delta;
  delta = sqrt(square(a - d) + 4 * b * c);

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B(2, 2);
  T half_delta = 0.5 * delta;
  T cosh_half_delta = cosh(half_delta);
  T sinh_half_delta = sinh(half_delta);
  T exp_half_a_plus_d = exp(0.5 * (a + d));
  T Two_exp_sinh = 2 * exp_half_a_plus_d * sinh_half_delta;
  T delta_cosh = delta * cosh_half_delta;
  T ad_sinh_half_delta = (a - d) * sinh_half_delta;

  B(0, 0) = exp_half_a_plus_d * (delta_cosh + ad_sinh_half_delta);
  B(0, 1) = b * Two_exp_sinh;
  B(1, 0) = c * Two_exp_sinh;
  B(1, 1) = exp_half_a_plus_d * (delta_cosh - ad_sinh_half_delta);

  return B / delta;
}

}  // namespace math
}  // namespace stan

#endif
