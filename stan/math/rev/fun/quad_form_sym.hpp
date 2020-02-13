#ifndef STAN_MATH_REV_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_REV_FUN_QUAD_FORM_SYM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/fun/quad_form.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * Symmetry of the resulting matrix is guaranteed.
 *
 * @tparam Ta type of elements in the symmetric matrix
 * @tparam Ra number of rows in the symmetric matrix, can be Eigen::Dynamic
 * @tparam Ca number of columns in the symmetric matrix, can be Eigen::Dynamic
 * @tparam Tb type of elements in the second matrix
 * @tparam Rb number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam Cb number of columns in the second matrix, can be Eigen::Dynamic
 *
 * @param A symmetric matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix of size Cb.
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb,
          require_any_var_t<Ta, Tb>...>
inline Eigen::Matrix<var, Cb, Cb> quad_form_sym(
    const Eigen::Matrix<Ta, Ra, Ca>& A, const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_symmetric("quad_form_sym", "A", A);
  check_multiplicable("quad_form_sym", "A", A, "B", B);

  internal::quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>* baseVari
      = new internal::quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>(A, B, true);

  return baseVari->impl_->C_;
}

/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * @tparam Ta type of elements in the symmetric matrix
 * @tparam Ra number of rows in the symmetric matrix, can be Eigen::Dynamic
 * @tparam Ca number of columns in the symmetric matrix, can be Eigen::Dynamic
 * @tparam Tb type of elements in the vector
 * @tparam Rb number of rows in the vector, can be Eigen::Dynamic
 *
 * @param A symmetric matrix
 * @param B vector
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename Ta, int Ra, int Ca, typename Tb, int Rb,
          require_any_var_t<Ta, Tb>...>
inline var quad_form_sym(const Eigen::Matrix<Ta, Ra, Ca>& A,
                         const Eigen::Matrix<Tb, Rb, 1>& B) {
  check_symmetric("quad_form_sym", "A", A);
  check_multiplicable("quad_form_sym", "A", A, "B", B);

  internal::quad_form_vari<Ta, Ra, Ca, Tb, Rb, 1>* baseVari
      = new internal::quad_form_vari<Ta, Ra, Ca, Tb, Rb, 1>(A, B, true);

  return baseVari->impl_->C_(0, 0);
}

}  // namespace math
}  // namespace stan
#endif
