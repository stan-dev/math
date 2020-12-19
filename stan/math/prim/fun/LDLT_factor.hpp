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
 * LDLT_factor is a thin wrapper on Eigen::LDLT to allow for
 * reusing factorizations and efficient autodiff of things like
 * log determinants and solutions to linear systems.
 *
 * Memory is allocated in the constructor and stored in a
 * <code>std::shared_ptr</code>, which ensures that is freed
 * when the object is released.
 *
 * After the constructor and/or compute() is called, users of
 * LDLT_factor are responsible for calling success() to
 * check whether the factorization has succeeded.  Use of an LDLT_factor
 * object (e.g., in mdivide_left_ldlt) is undefined if success() is false.
 *
 * Its usage pattern is:
 *
 * ~~~
 * Eigen::Matrix<T, R, C> A1, A2;
 *
 * LDLT_factor<T, R, C> ldlt_A1(A1);
 * LDLT_factor<T, R, C> ldlt_A2;
 * ldlt_A2.compute(A2);
 * ~~~
 *
 * The caller should check that ldlt_A1.success() and ldlt_A2.success()
 * are true or abort accordingly.  Alternatively, call check_ldlt_factor().
 *
 * Note that ldlt_A1 and ldlt_A2 are completely equivalent.  They simply
 * demonstrate two different ways to construct the factorization.
 *
 * The caller can use the LDLT_factor objects as needed.  For
 * instance
 *
 * ~~~
 * x1 = mdivide_left_ldlt(ldlt_A1, b1);
 * x2 = mdivide_right_ldlt(b2, ldlt_A2);
 *
 * d1 = log_determinant_ldlt(ldlt_A1);
 * d2 = log_determinant_ldlt(ldlt_A2);
 * ~~~
 *
 * This class is conceptually similar to the corresponding Eigen
 * class.  Any symmetric, positive-definite matrix A can be
 * decomposed as LDL' where L is unit lower-triangular and D is
 * diagonal with positive diagonal elements.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 */
template <typename T, bool alloc_in_arena>
class LDLT_factor2;

template <typename T>
class LDLT_factor2<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, false> {
private:
  using ldlt_type = Eigen::LDLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
  std::shared_ptr<ldlt_type> ldlt_ptr_;
public:
  template <typename S,
	    require_eigen_t<S>* = nullptr>
  LDLT_factor2(const S& matrix) :
    ldlt_ptr_(new ldlt_type(matrix)) {}

  template <typename Rhs>
  inline const Eigen::Solve<ldlt_type, Rhs> solve(
      const Eigen::MatrixBase<Rhs>& b) const {
    return ldlt_ptr_->solve(b);
  }

  auto& matrix() const {
    return ldlt_ptr_->matrixLDLT();
  }

  auto& ldlt() const {
    return *ldlt_ptr_;
  }
};

template <typename T, typename... Args,
	  require_matrix_t<T>* = nullptr>
inline auto make_ldlt_factor(const T& A) {
  return LDLT_factor2<plain_type_t<T>,
		      disjunction<is_var<scalar_type_t<return_type_t<T, Args...>>>,
				  conjunction<is_fvar<return_type_t<Args...>>,
					      is_var<partials_type_t<return_type_t<Args...>>>>,
				  conjunction<is_fvar<return_type_t<Args...>>,
					      is_fvar<partials_type_t<return_type_t<Args...>>>,
					      is_var<partials_type_t<partials_type_t<return_type_t<Args...>>>>>>::value>(A);
}

  
}  // namespace math
}  // namespace stan

#endif
