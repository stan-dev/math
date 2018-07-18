#ifndef STAN_MATH_PRIM_MAT_FUN_LDLT_FACTOR_HPP
#define STAN_MATH_PRIM_MAT_FUN_LDLT_FACTOR_HPP

#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>

namespace stan {
namespace math {

/**
 * LDLT_factor is a wrapper on Eigen::LDLT which factorizes the matrix A into P
 * * L * D * L^T * P^T where P is a permutation matrix (permuting the rows of
 * the matrix it is multiplied with), L is a lower triangular matrix and, D is a
 * diagonal matrix
 *
 * In this implementation, the matrix A must be positive definite
 *
 * The two usage patterns are:
 *
 * Eigen::Matrix<T, R, C> A;
 * LDLT_factor<T, R, C> ldlt1(A);
 *
 * which computes the factorization in the constructor and
 *
 * LDLT_factor<T, R, C> ldlt2;
 * ldlt2.compute(A);
 *
 * which defers the factorization until later.
 *
 * @tparam T scalar type of matrix
 * @tparam R row type of matrix (as in Eigen)
 * @tparam C column type of matrix (as in Eigen)
 */
template <typename T, int R, int C>
class LDLT_factor {
 private:
  int N_;
  Eigen::LDLT<Eigen::Matrix<T, R, C> > ldlt_;
  bool success_;

 public:
  LDLT_factor() : N_(0), success_(false) {}

  /**
   * Construct a new LDLT_factor that contains the LDLT factorization of the
   * matrix A
   *
   * @param A matrix to factorize
   */
  explicit LDLT_factor(const Eigen::Matrix<T, R, C>& A)
      : N_(A.rows()), success_(false) {
    compute(A);
  }

  /**
   * @return number of rows of factorized matrix
   */
  inline int rows() const { return N_; }

  /**
   * @return number of columns of factorized matrix
   */
  inline int cols() const { return N_; }

  /**
   * Compute the LDLT factorization of the matrix A
   *
   * @param A matrix to factorize
   * @throw std::domain_error if the Eigen LDLT factorization fails, A is not
   * positive-semidefinite, or the diagonal in the factorization contains NaNs
   */
  inline void compute(const Eigen::Matrix<T, R, C>& A) {
    check_square("LDLT_factor", "A", A);

    N_ = A.rows();
    ldlt_.compute(A);
    success_ = false;

    if (ldlt_.info() != Eigen::Success) {
      domain_error("compute", "A", "", "Eigen failed to compute ldlt for A",
                   "");
    }

    if (!ldlt_.isPositive()) {
      domain_error("compute", "A", "", "is not positive semidefinite", "");
    }

    Eigen::Matrix<T, Eigen::Dynamic, 1> ldlt_diag(ldlt_.vectorD());
    for (int i = 0; i < ldlt_diag.size(); ++i) {
      if (is_nan(ldlt_diag(i))) {
        domain_error("compute", "Diagonal of factorization", "", "contains NaN",
                     "");
      }

      if (ldlt_diag(i) <= 0.0) {
        domain_error("compute", "Diagonal of factorization", "",
                     "should be positive", "");
      }
    }

    success_ = true;
  }

  /**
   * Return the log of the absolute value of the determinant
   * of the factorized matrix
   *
   * @return log(abs(determinant(A)))
   */
  inline T log_abs_det() const { return ldlt_.vectorD().array().log().sum(); }

  /**
   * Solve the system A x = b for x
   *
   * Store the result in dst. dst should already be the right
   * size to hold the solution
   *
   * @tparam RhsType type of right hand side
   * @tparam DstType type of solution
   * @param rhs right hand side of equation (b)
   * @param dst solution (x)
   */
  template <typename RhsType, typename DstType>
  inline void solve(const Eigen::MatrixBase<RhsType>& rhs,
                    Eigen::MatrixBase<DstType>& dst) const {
    dst = ldlt_.solve(rhs);
  }

  /**
   * Solve the system A x = b for x and return the result
   *
   * @tparam RhsType type of right hand side
   * @param rhs right hand side of equation (b)
   * @return solution (x)
   */
  template <typename RhsType>
  inline Eigen::Matrix<T, R, C> solve(
      const Eigen::MatrixBase<RhsType>& rhs) const {
    return ldlt_.solve(rhs);
  }

  /**
   * Solve the system A x = b for x in place
   *
   * @tparam RhsType Type of right hand side
   * @param[in, out] rhs on input rhs is b and on output it is x
   */
  template <typename RhsType>
  inline void solveInPlace(Eigen::MatrixBase<RhsType>& rhs) const {
    ldlt_.solveInPlace(rhs);
  }

  /*
   * Return true if:
   *  1. Eigen reported no errors in factorization
   *  2. Matrix is positive semi-definite
   *  3. There are no NaNs in the diagnonal of the factorization
   *
   * @return whether or not decomposition was a success
   */
  inline bool success() const { return success_; }

  /*
   * Return the non-zero values of diagonal matrix D in the LDL^T decomposition
   * of A
   *
   * @return non-zero values of D
   */
  inline Eigen::Matrix<T, Eigen::Dynamic, 1> vectorD() const {
    return ldlt_.vectorD();
  }
};

}  // namespace math
}  // namespace stan
#endif
