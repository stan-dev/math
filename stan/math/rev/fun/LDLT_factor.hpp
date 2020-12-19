#ifndef STAN_MATH_REV_FUN_LDLT_FACTOR_HPP
#define STAN_MATH_REV_FUN_LDLT_FACTOR_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>

namespace stan {
namespace math {

/**
 * An LDLT_factor with `alloc_in_arena = True` is a structure that
 * holds a matrix of type T and the LDLT of its values, where all
 * member variable allocations are done in the arena.
 *
 * @tparam T type of elements in the matrix
 */
template <typename T>
class LDLT_factor<T, true> {
private:
  using ldlt_type = Eigen::LDLT<T>;
  ldlt_type* ldlt_ptr_;
public:
  template <typename S,
	    require_eigen_t<S>* = nullptr>
  LDLT_factor(const S& matrix) :
    ldlt_ptr_(make_chainable_ptr(matrix.ldlt())) {}

  /**
   * Return a const reference to the underlying matrix
   */
  const auto& matrix() const {
    return ldlt_ptr_->matrixLDLT();
  }

  /**
   * Return a const reference to the LDLT factor of the matrix
   */
  const auto& ldlt() const {
    return *ldlt_ptr_;
  }
};
  
/**
 * An LDLT_factor of an `Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>`
 * with `alloc_in_arena = True` holds a copy of the input matrix and the LDLT
 * of its values, with all member variable allocations are done in the arena.
 */
template <>
class LDLT_factor<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>, true> {
private:
  arena_t<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> matrix_;
  Eigen::LDLT<Eigen::MatrixXd>* ldlt_ptr_;
public:
  template <typename S,
	    require_eigen_vt<is_var, S>* = nullptr>
  LDLT_factor(const S& matrix) :
    matrix_(matrix),
    ldlt_ptr_(make_chainable_ptr(matrix.val().ldlt())) {}

  /**
   * Return the underlying matrix
   */
  auto matrix() const {
    return matrix_;
  }

  /**
   * Return a const reference to the LDLT factor of the matrix
   */
  const auto& ldlt() const {
    return *ldlt_ptr_;
  }
};

/**
 * An LDLT_factor of a `var_value<Eigen::MatrixXd>` with `alloc_in_arena = True`
 * holds a copy of the input `var_value` and the LDLT of its values,
 * with all member variable allocations are done in the arena.
 */
template <>
class LDLT_factor<var_value<Eigen::MatrixXd>, true> {
private:
  var_value<Eigen::MatrixXd> matrix_;
  Eigen::LDLT<Eigen::MatrixXd>* ldlt_ptr_;
public:
  template <typename S,
	    require_var_matrix_t<S>* = nullptr>
  LDLT_factor(const S& matrix) :
    matrix_(matrix), ldlt_ptr_(make_chainable_ptr(matrix.val().ldlt())) {}

  /**
   * Return a const reference the underlying `var_value`
   */
  const auto& matrix() const {
    return matrix_;
  }

  /**
   * Return a const reference to the LDLT factor of the matrix
   */
  const auto& ldlt() const {
    return *ldlt_ptr_;
  }
};

}  // namespace math
}  // namespace stan
#endif
