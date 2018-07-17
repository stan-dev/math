#ifndef STAN_MATH_REV_MAT_FUN_LDLT_FACTOR_HPP
#define STAN_MATH_REV_MAT_FUN_LDLT_FACTOR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/LDLT_factor.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <limits>

namespace stan {
namespace math {

/**
 * LDLT_factor is a wrapper on Eigen::LDLT which factorizes the matrix A into P
 * * L * D * L^T * P^T where P is a matrix that reorders the rows of L (for
 * numerical stability), L is a lower triangular matrix and, D is a diagonal
 * matrix
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
template <int R, int C>
class LDLT_factor<var, R, C> {
 private:
  int N_;
  int *transpositions_mem_;  // this memory holds P
  double *vectorD_mem_;
  double *ldlt_mem_;
  vari **variA_mem_;
  bool success_;

 public:
  LDLT_factor()
      : N_(0),
        transpositions_mem_(NULL),
        vectorD_mem_(NULL),
        ldlt_mem_(NULL),
        variA_mem_(NULL),
        success_(false) {}

  /**
   * Construct a new LDLT_factor that contains the LDLT factorization of the
   * matrix A
   *
   * @param A matrix to factorize
   */
  explicit LDLT_factor(const Eigen::Matrix<var, R, C> &A) : N_(A.rows()) {
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
   * Use the LDLT_factor object to factorize a new matrix.  After calling
   * this function, the user should call success() to check that the
   * factorization was successful. If the factorization is not successful,
   * the LDLT_factor is not valid and other functions should not be used.
   *
   * @param A A symmetric positive definite matrix to factorize
   */
  inline void compute(const Eigen::Matrix<var, R, C> &A) {
    check_square("compute", "A", A);

    ldlt_mem_
        = ChainableStack::instance().memalloc_.alloc_array<double>(A.size());
    variA_mem_
        = ChainableStack::instance().memalloc_.alloc_array<vari *>(A.size());
    Eigen::Map<Eigen::MatrixXd> ldlt_matrix(ldlt_mem_, A.rows(), A.cols());

    N_ = A.rows();
    success_ = false;

    for (int j = 0; j < A.cols(); j++) {
      for (int i = 0; i < A.rows(); i++) {
        ldlt_matrix(i, j) = A(i, j).val();
        variA_mem_[i + N_ * j] = A(i, j).vi_;
      }
    }

    Eigen::LDLT<Eigen::Ref<Eigen::Matrix<double, R, C> > > ldlt(ldlt_matrix);

    const auto &transpositions_indices = ldlt.transpositionsP().indices();
    transpositions_mem_
        = ChainableStack::instance().memalloc_.alloc_array<int>(N_);
    for (int i = 0; i < N_; ++i) {
      transpositions_mem_[i] = transpositions_indices(i);
    }

    const auto &vectorD = ldlt.vectorD();
    vectorD_mem_ = ChainableStack::instance().memalloc_.alloc_array<double>(N_);
    for (int i = 0; i < N_; ++i) {
      vectorD_mem_[i] = vectorD(i);
    }

    if (ldlt.info() != Eigen::Success) {
      domain_error("compute", "A", "", "Eigen failed to compute ldlt for A",
                   "");
    }

    if (!ldlt.isPositive()) {
      domain_error("compute", "A", "", "is not positive semidefinite", "");
    }

    Eigen::VectorXd ldlt_diag(ldlt.vectorD());
    for (int i = 0; i < ldlt_diag.size(); ++i) {
      if (is_nan(ldlt_diag(i))) {
        domain_error("compute", "A", "", "contains NaN in factorization", "");
      }

      if (ldlt_diag(i) <= 0.0) {
        domain_error("compute", "Diagonal of factorization", "",
                     "should be positive", "");
      }
    }

    success_ = true;
  }

  /**
   * Access the varis from the input matrix A
   *
   * @param i row of vari
   * @param j column of vari
   * @return vari
   */
  vari *getVariA(int i, int j) { return variA_mem_[i + N_ * j]; }

  /**
   * Compute the log of the absolute value of the determinant
   * of the factorized matrix
   *
   * @return log(abs(determinant(A)))
   */
  inline double log_abs_det() const { return vectorD().array().log().sum(); }

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
  void solve(const Eigen::MatrixBase<RhsType> &rhs,
             Eigen::MatrixBase<DstType> &dst) const {
    check_size_match("solve", "Rows of ", "right hand side", rhs.rows(),
                     "rows of ", "solution vector", dst.rows());

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
        ldlt_matrix(ldlt_mem_, N_, N_);
    Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, 1> > transposition_indices(
        transpositions_mem_, N_);
    Eigen::Transpositions<Eigen::Dynamic> m_transpositions(
        transposition_indices);

    dst = m_transpositions * rhs;
    ldlt_matrix.triangularView<Eigen::UnitLower>().solveInPlace(dst);

    for (int i = 0; i < N_; ++i) {
      if (std::abs(vectorD_mem_[i]) > std::numeric_limits<double>::min())
        dst.row(i) /= vectorD_mem_[i];
      else
        dst.row(i).setZero();
    }

    ldlt_matrix.triangularView<Eigen::UnitLower>().transpose().solveInPlace(
        dst);
    dst = m_transpositions.transpose() * dst;
  }

  /**
   * Solve the system A x = b for x and return the result
   *
   * @tparam RhsType type of right hand side
   * @param rhs right hand side of equation (b)
   * @return solution (x)
   */
  template <typename RhsType>
  inline Eigen::Matrix<double, R, C> solve(
      const Eigen::MatrixBase<RhsType> &rhs) const {
    Eigen::Matrix<double, R, C> dst(rhs.rows(), rhs.cols());
    solve(rhs, dst);
    return dst;
  }

  /**
   * Solve the system A x = b for x in place
   *
   * @tparam RhsType Type of right hand side
   * @param[in, out] rhs on input rhs is b and on output it is x
   */
  template <typename RhsType>
  inline void solveInPlace(Eigen::MatrixBase<RhsType> &rhs) const {
    solve(rhs, rhs);
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
  inline Eigen::VectorXd vectorD() const {
    Eigen::VectorXd vectorD(N_);
    for (int i = 0; i < N_; ++i)
      vectorD(i) = vectorD_mem_[i];
    return vectorD;
  }
};

}  // namespace math
}  // namespace stan
#endif
