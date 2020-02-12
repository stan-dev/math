#ifndef STAN_MATH_FWD_FUN_QUAD_FORM_HPP
#define STAN_MATH_FWD_FUN_QUAD_FORM_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/dot_product.hpp>
#include <stan/math/fwd/fun/multiply.hpp>

namespace stan {
namespace math {

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * @tparam Ta type of elements in the square matrix
 * @tparam Ra number of rows in the square matrix, can be Eigen::Dynamic
 * @tparam Ca number of columns in the square matrix, can be Eigen::Dynamic
 * @tparam Tb type of elements in the second matrix
 * @tparam Rb number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam Cb number of columns in the second matrix, can be Eigen::Dynamic
 *
 * @param A square matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix of size Cb.
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb,
          require_any_fvar_t<Ta, Tb>...>
inline Eigen::Matrix<return_type_t<Ta, Tb>, Cb, Cb> quad_form(
    const Eigen::Matrix<Ta, Ra, Ca>& A, const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);
  return multiply(transpose(B), multiply(A, B));
}

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * @tparam Ta type of elements in the square matrix
 * @tparam Ra number of rows in the square matrix, can be Eigen::Dynamic
 * @tparam Ca number of columns in the square matrix, can be Eigen::Dynamic
 * @tparam Tb type of elements in the vector
 * @tparam Rb number of rows in the vector, can be Eigen::Dynamic
 *
 * @param A square matrix
 * @param B vector
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename Ta, int Ra, int Ca, typename Tb, int Rb,
          require_any_fvar_t<Ta, Tb>...>
inline return_type_t<Ta, Tb> quad_form(const Eigen::Matrix<Ta, Ra, Ca>& A,
                                       const Eigen::Matrix<Tb, Rb, 1>& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);
  return dot_product(B, multiply(A, B));
}

}  // namespace math
}  // namespace stan

#endif
