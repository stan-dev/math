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
template <typename T, int R, int C>
class LDLT_factor {
 public:
  using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  using matrix_t = Eigen::Matrix<T, R, C>;
  using ldlt_t = Eigen::LDLT<matrix_t>;
  using size_type = size_t;
  using value_type = double;

  LDLT_factor() : N_(0), ldltP_(new ldlt_t()) {}

  LDLT_factor(const matrix_t& A) : N_(0), ldltP_(new ldlt_t()) { compute(A); }

  inline void compute(const matrix_t& A) {
    check_square("LDLT_factor", "A", A);
    N_ = A.rows();
    ldltP_->compute(A);
  }

  inline bool success() const {
    if (ldltP_->info() != Eigen::Success) {
      return false;
    }
    if (!(ldltP_->isPositive())) {
      return false;
    }
    vector_t ldltP_diag(ldltP_->vectorD());
    for (int i = 0; i < ldltP_diag.size(); ++i) {
      if (ldltP_diag(i) <= 0 || is_nan(ldltP_diag(i))) {
        return false;
      }
    }
    return true;
  }

  inline T log_abs_det() const { return sum(log(ldltP_->vectorD())); }

  inline void inverse(matrix_t& invA) const {
    invA.setIdentity(N_);
    ldltP_->solveInPlace(invA);
  }

  template <typename Rhs>
  inline const Eigen::Solve<ldlt_t, Rhs> solve(
      const Eigen::MatrixBase<Rhs>& b) const {
    return ldltP_->solve(b);
  }

  inline matrix_t solveRight(const matrix_t& B) const {
    return ldltP_->solve(B.transpose()).transpose();
  }

  inline vector_t vectorD() const { return ldltP_->vectorD(); }

  // This could return an ldlt_t, but that would just do another ldlt on a
  // matrix constructed from an ldlt
  inline matrix_t matrixLDLT() const { return ldltP_->matrixLDLT(); }

  inline size_t rows() const { return N_; }
  inline size_t cols() const { return N_; }

  size_t N_;
  std::shared_ptr<ldlt_t> ldltP_;
};

template <typename T, bool alloc_in_arena>
class LDLT_factor2;

template <typename T>
class LDLT_factor2<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, false> {
 private:
  using ldlt_type
      = Eigen::LDLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
  std::shared_ptr<ldlt_type> ldlt_ptr_;

 public:
  template <typename S, require_eigen_t<S>* = nullptr>
  LDLT_factor2(const S& matrix) : ldlt_ptr_(new ldlt_type(matrix)) {}

  template <typename Rhs>
  inline const Eigen::Solve<ldlt_type, Rhs> solve(
      const Eigen::MatrixBase<Rhs>& b) const {
    return ldlt_ptr_->solve(b);
  }

  auto& matrix() { return ldlt_ptr_->matrixLDLT(); }

  auto& ldlt() { return *ldlt_ptr_; }
};

/*template <typename T,
          // require_not_vt_var<T>* = nullptr, <-- using this instead of next
line doesn't work, why?
          //require_not_eigen_vt<is_var, T>* = nullptr,
        require_not_rev_matrix_t<T>* = nullptr>
inline auto make_ldlt_factor(const T& A) {
return LDLT_factor<value_type_t<T>, T::RowsAtCompileTime,
T::ColsAtCompileTime>(A);
}*/

template <typename T, typename... Args, require_matrix_t<T>* = nullptr>
inline auto make_ldlt_factor(const T& A) {
  return LDLT_factor2<
      plain_type_t<T>,
      disjunction<
          is_var<scalar_type_t<return_type_t<T, Args...>>>,
          std::is_same<fvar<var>, return_type_t<Args...>>,
          std::is_same<fvar<fvar<var>>, return_type_t<Args...>>>::value>(A);
}

}  // namespace math
}  // namespace stan

#endif
