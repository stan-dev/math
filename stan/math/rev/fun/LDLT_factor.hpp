#ifndef STAN_MATH_REV_FUN_LDLT_FACTOR_HPP
#define STAN_MATH_REV_FUN_LDLT_FACTOR_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_alloc.hpp>
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
class LDLT_factor<T, require_rev_matrix_t<T>> {
  arena_t<T> A_;
  decltype(make_chainable_ptr(std::declval<T>().val().ldlt())) ldlt_ptr_;
public:
  explicit LDLT_factor(const T &A)
    : A_(A), ldlt_ptr_(make_chainable_ptr(A.val().ldlt())) {
  }

  template <typename S>
  inline auto solve(const S& b) {
    return ldlt_ptr_->get().solve(b).eval();
  }

  template <typename S>
  inline void solveInPlace(S& b) {
    return ldlt_ptr_->get().solveInPlace(b);
  }

  inline bool success() const {
    return A_.rows() != 0 &&
      ldlt_ptr_->get().info() == Eigen::Success &&
      ldlt_ptr_->get().isPositive() &&
      (ldlt_ptr_->get().vectorD().array() > 0).all();
  }

  inline double log_abs_det() const { return sum(log(vectorD())); }

  inline Eigen::VectorXd vectorD() {
    return ldlt_ptr_->get().vectorD();
  }

  inline size_t rows() const {
    return A_.rows();
  }

  inline size_t cols() const {
    return A_.cols();
  }

  using size_type = size_t;
  using value_type = var;
};

}  // namespace math
}  // namespace stan
#endif
