#ifndef STAN_MATH_REV_FUN_LDLT_FACTOR_HPP
#define STAN_MATH_REV_FUN_LDLT_FACTOR_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>

namespace stan {
namespace math {

/**
 * An LDLT_factor of an `Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>`
 * with `alloc_in_arena = True` holds a copy of the input matrix and the LDLT
 * of its values, with all member variable allocations are done in the arena.
 */
template <typename T>
class LDLT_factor<T,
		  std::enable_if_t<bool_constant<is_eigen_matrix_dynamic<T>::value &&
						 is_var<scalar_type_t<T>>::value>::value>> {
 private:
  const plain_type_t<T>& matrix_;
  Eigen::LDLT<Eigen::MatrixXd> ldlt_;

 public:
  template <typename S, require_same_t<plain_type_t<T>, plain_type_t<S>>* = nullptr>
  explicit LDLT_factor(const S& matrix)
      : matrix_(matrix), ldlt_(matrix.val().ldlt()) {}

  /**
   * Return a const reference to the underlying matrix
   */
  const auto& matrix() const { return matrix_; }

  /**
   * Return a const reference to the LDLT factor of the matrix values
   */
  const auto& ldlt() const { return ldlt_; }
};

/**
 * An LDLT_factor of a `var_value<Eigen::MatrixXd>`
 * holds a copy of the input `var_value` and the LDLT of its values.
 */
template <typename T>
class LDLT_factor<T,
		  std::enable_if_t<is_var_matrix<T>::value>> {
 private:
  plain_type_t<T> matrix_;
  Eigen::LDLT<Eigen::MatrixXd> ldlt_;

 public:
  template <typename S, require_same_t<plain_type_t<T>, plain_type_t<S>>* = nullptr>
  explicit LDLT_factor(const S& matrix)
      : matrix_(matrix), ldlt_(matrix.val().ldlt()) {}

  /**
   * Return a const reference the underlying `var_value`
   */
  const auto& matrix() const { return matrix_; }

  /**
   * Return a const reference to the LDLT factor of the matrix values
   */
  const auto& ldlt() const { return ldlt_; }
};

}  // namespace math
}  // namespace stan
#endif
