#ifndef STAN_MATH_REV_FUN_LDLT_FACTOR_HPP
#define STAN_MATH_REV_FUN_LDLT_FACTOR_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>

namespace stan {
namespace math {

/**
 * A template specialization of src/stan/math/matrix/LDLT_factor.hpp for
 * var which can be used with all the *_ldlt functions.
 *
 * The usage pattern is:
 *
 * ~~~
 * Eigen::Matrix<T, R, C> A1, A2;
 *
 * LDLT_factor<T, R, C> ldlt_A1(A1);
 * LDLT_factor<T, R, C> ldlt_A2;
 * ldlt_A2.compute(A2);
 * ~~~
 *
 * Now, the caller should check that ldlt_A1.success() and ldlt_A2.success()
 * are true or abort accordingly.  Alternatively, call check_ldlt_factor().
 * The behaviour of using an LDLT_factor without success() returning true is
 * undefined.
 *
 * Note that ldlt_A1 and ldlt_A2 are completely equivalent.  They simply
 * demonstrate two different ways to construct the factorization.
 *
 * Now, the caller can use the LDLT_factor objects as needed.  For instance
 *
 * ~~~
 * x1 = mdivide_left_ldlt(ldlt_A1, b1);
 * x2 = mdivide_right_ldlt(b2, ldlt_A2);
 *
 * d1 = log_determinant_ldlt(ldlt_A1);
 * d2 = log_determinant_ldlt(ldlt_A2);
 * ~~~
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 */

template <typename T>
class LDLT_factor2<T, true> {
private:
  using ldlt_type = Eigen::LDLT<T>;
  ldlt_type* ldlt_ptr_;
public:
  template <typename S,
	    require_eigen_t<S>* = nullptr>
  LDLT_factor2(const S& matrix) :
    ldlt_ptr_(make_chainable_ptr(matrix.ldlt())) {}

  const auto& matrix() const {
    return ldlt_ptr_->matrixLDLT();
  }

  const auto& ldlt() const {
    return *ldlt_ptr_;
  }
};
  
template <>
class LDLT_factor2<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>, true> {
private:
  arena_t<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> matrix_;
  Eigen::LDLT<Eigen::MatrixXd>* ldlt_ptr_;
public:
  template <typename S,
	    require_eigen_vt<is_var, S>* = nullptr>
  LDLT_factor2(const S& matrix) :
    matrix_(matrix),
    ldlt_ptr_(make_chainable_ptr(matrix.val().ldlt())) {}

  auto matrix() const {
    return matrix_;
  }

  const auto& ldlt() const {
    return *ldlt_ptr_;
  }
};

template <>
class LDLT_factor2<var_value<Eigen::MatrixXd>, true> {
private:
  var_value<Eigen::MatrixXd> matrix_;
  Eigen::LDLT<Eigen::MatrixXd>* ldlt_ptr_;
public:
  template <typename S,
	    require_var_matrix_t<S>* = nullptr>
  LDLT_factor2(const S& matrix) :
    matrix_(matrix), ldlt_ptr_(make_chainable_ptr(matrix.val().ldlt())) {}

  auto matrix() const {
    return matrix_;
  }

  const auto& ldlt() const {
    return *ldlt_ptr_;
  }
};

}  // namespace math
}  // namespace stan
#endif
