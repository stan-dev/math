#ifndef STAN_MATH_PRIM_FUN_LDLT_FACTOR_HPP
#define STAN_MATH_PRIM_FUN_LDLT_FACTOR_HPP

#include <stan/math/prim/err/check_square.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <memory>

namespace stan {
namespace math {

/**
 * LDLT_factor is a structure that holds a matrix of type T and the
 * LDLT of its values.
 *
 * @tparam T type of elements in the matrix
 */
template <typename T, typename Enable = void>
class LDLT_factor;

/**
 * An LDLT_factor is a structure that holds a matrix of type T and the
 * LDLT of its values.
 *
 * @tparam T type of elements in the matrix
 */
template <typename T>
class LDLT_factor<T, std::enable_if_t<bool_constant<
                         is_eigen_matrix_dynamic<T>::value
                         && !is_var<scalar_type_t<T>>::value>::value>> {
 private:
  plain_type_t<T> matrix_;
  Eigen::LDLT<plain_type_t<T>> ldlt_;

 public:
  template <typename S,
            require_same_t<plain_type_t<T>, plain_type_t<S>>* = nullptr>
  explicit LDLT_factor(const S& matrix)
      : matrix_(matrix), ldlt_(matrix_.ldlt()) {}

  /**
   * Return a const reference to the underlying matrix
   */
  const auto& matrix() const { return matrix_; }

  /**
   * Return a const reference to the LDLT factor of the matrix
   */
  const auto& ldlt() const { return ldlt_; }
};

/**
 * Make an LDLT_factor with matrix type `plain_type_t<T>`
 *
 * @tparam T Type of matrix to take the LDLT of
 * @param A Matrix to take the LDLT of
 * @return LDLT_factor of A
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline auto make_ldlt_factor(const T& A) {
  return LDLT_factor<plain_type_t<T>>(A);
}

}  // namespace math
}  // namespace stan

#endif
