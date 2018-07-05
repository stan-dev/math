#ifndef STAN_MATH_REV_MAT_FUN_LDLT_FACTOR_HPP
#define STAN_MATH_REV_MAT_FUN_LDLT_FACTOR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/LDLT_factor.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

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
 */
template <int R, int C>
class LDLT_factor<var, R, C> {
 private:
  int N_;
  Eigen::LDLT<Eigen::Ref<Eigen::Map<Eigen::Matrix<double, R, C> > > > *ldltP_;
  vari **variA_mem_;

 public:
  /**
   * Default constructor.  The caller *MUST* call compute() after this.  Any
   * calls which use the LDLT_factor without calling compute() run the risk
   * of crashing Stan from within Eigen.
   */
  LDLT_factor() : N_(0), ldltP_(NULL), variA_mem_(NULL) {}

  explicit LDLT_factor(const Eigen::Matrix<var, R, C> &A) { compute(A); }

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

    void *ldlt_class_mem = ChainableStack::instance().memalloc_.alloc(sizeof(
        Eigen::LDLT<Eigen::Ref<Eigen::Map<Eigen::Matrix<double, R, C> > > >));
    double *ldlt_mem
        = ChainableStack::instance().memalloc_.alloc_array<double>(A.size());
    variA_mem_
        = ChainableStack::instance().memalloc_.alloc_array<vari *>(A.size());
    Eigen::Map<Eigen::MatrixXd> ldlt_matrix(ldlt_mem, A.rows(), A.cols());

    N_ = A.rows();

    for (int j = 0; j < A.cols(); j++) {
      for (int i = 0; i < A.rows(); i++) {
        ldlt_matrix(i, j) = A(i, j).val();
        variA_mem_[i + N_ * j] = A(i, j).vi_;
      }
    }

    ldltP_ = new (ldlt_class_mem)
        Eigen::LDLT<Eigen::Ref<Eigen::Map<Eigen::Matrix<double, R, C> > > >(
            ldlt_matrix);
  }

  /**
   * Access the varis from the input matrix A
   *
   * @param i row of vari
   * @param j column of vari
   * @return vari
   */
  vari *get_variA(int i, int j) { return variA_mem_[i + N_ * j]; }

  /**
   * Compute the log of the absolute value of the determinant
   * of the input matrix A
   *
   * @return log(abs(det(A)))
   */
  inline double log_abs_det() const { return vectorD().array().log().sum(); }

  /**
   * Compute the actual numerical result of inv(A)*b.  Note that this isn't
   * meant to handle any of the autodiff.  This is a convenience function
   * for the actual implementations in mdivide_left_ldlt.
   *
   * Precondition: success() must return true. If success() returns false,
   *    this function runs the risk of crashing Stan from within Eigen.
   *
   * @param b The right handside.  Note that this is templated such that
   * Eigen's expression-templating magic can work properly here.
   */
  template <typename Rhs>
  inline const Eigen::Solve<
      Eigen::LDLT<Eigen::Ref<Eigen::Map<Eigen::Matrix<double, R, C> > > >, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    return ldltP_->solve(b);
  }

  /**
   * Solve the system inv(A) * b in place (store solution where b was
   * originally)
   *
   * template Rhs Type of right hand side
   * @param b Right hand side of equation
   */
  template <typename Rhs>
  inline void solveInPlace(Eigen::MatrixBase<Rhs> &b) const {
    ldltP_->solveInPlace(b);
  }

  /**
   * Determine whether the most recent factorization succeeded.  This should
   * always be called after the object is constructed (with a matrix) or
   * after compute() is called.
   */
  inline bool success() const {
    bool ret;
    ret = N_ != 0;
    ret = ret && ldltP_->info() == Eigen::Success;
    ret = ret && ldltP_->isPositive();
    ret = ret && (ldltP_->vectorD().array() > 0).all();
    return ret;
  }

  /**
   * The entries of the diagonal matrix D.  They should be strictly positive
   * for a positive definite matrix.
   *
   * Precondition: success() must return true. If success() returns false,
   *    this function runs the risk of crashing Stan from within Eigen.
   */
  inline Eigen::VectorXd vectorD() const { return ldltP_->vectorD(); }

  inline int rows() const { return N_; }
  inline int cols() const { return N_; }
};

}  // namespace math
}  // namespace stan
#endif
