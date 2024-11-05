#ifndef STAN_MATH_PRIM_FUN_EIGEN_HPP
#define STAN_MATH_PRIM_FUN_EIGEN_HPP

#ifdef EIGEN_MATRIXBASE_PLUGIN
#ifndef EIGEN_STAN_MATRIXBASE_PLUGIN
#error "Stan uses Eigen's EIGEN_MATRIXBASE_PLUGIN macro. To use your own "
"plugin add the eigen_plugin.h file to your plugin."
#endif
#else
#define EIGEN_MATRIXBASE_PLUGIN "stan/math/prim/eigen_plugins.h"
#endif

#ifdef EIGEN_ARRAYBASE_PLUGIN
#ifndef EIGEN_STAN_ARRAYBASE_PLUGIN
#error "Stan uses Eigen's EIGEN_ARRAYBASE_PLUGIN macro. To use your own "
    "plugin add the eigen_plugin.h file to your plugin."
#endif
#else
#define EIGEN_ARRAYBASE_PLUGIN "stan/math/prim/eigen_plugins.h"
#endif

namespace stan {
namespace math {
template <typename T>
class complex;
}
}
namespace Eigen {
namespace numext {
template <typename T>
inline auto& real(stan::math::complex<T>& x) {
  return x.real();
}

template <typename T>
inline auto& imag(stan::math::complex<T>& x) {
  return x.imag();
}

template <typename T>
inline const auto& real(const stan::math::complex<T>& x) {
  return x.real();
}

template <typename T>
inline const auto& imag(const stan::math::complex<T>& x) {
  return x.imag();
}
}

}
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/src/Core/NumTraits.h>
#include <Eigen/SVD>

namespace Eigen {
template <typename Real_>
struct NumTraits<stan::math::complex<Real_> > : GenericNumTraits<stan::math::complex<Real_> > {
  typedef Real_ Real;
  typedef typename NumTraits<Real_>::Literal Literal;
  enum {
    IsComplex = 1,
    IsInteger = 0,
    IsSigned = NumTraits<Real_>::IsSigned,
    RequireInitialization = NumTraits<Real_>::RequireInitialization,
    ReadCost = 2 * NumTraits<Real_>::ReadCost,
    AddCost = 2 * NumTraits<Real>::AddCost,
    MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
  };

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static inline Real epsilon() { return NumTraits<Real>::epsilon(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static inline Real dummy_precision() { return NumTraits<Real>::dummy_precision(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static inline int digits10() { return NumTraits<Real>::digits10(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static inline int max_digits10() { return NumTraits<Real>::max_digits10(); }
};
}
namespace Eigen {

  /**
   * Traits specialization for Eigen binary operations for `int`
   * and `double` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<int, double, BinaryOp> {
    using ReturnType = double;
  };

  /**
   * Traits specialization for Eigen binary operations for `double`
   * and `int` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<double, int, BinaryOp> {
    using ReturnType = double;
  };

  /**
   * Traits specialization for Eigen binary operations for `int`
   * and complex `double` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename T, typename BinaryOp>
  struct ScalarBinaryOpTraits<int, stan::math::complex<T>, BinaryOp> {
    using ReturnType = stan::math::complex<T>;
  };

  /**
   * Traits specialization for Eigen binary operations for complex
   * `double` and `int` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename T, typename BinaryOp>
  struct ScalarBinaryOpTraits<stan::math::complex<T>, int, BinaryOp> {
    using ReturnType = stan::math::complex<T>;
  };

  template <typename T, typename BinaryOp>
  struct ScalarBinaryOpTraits<stan::math::complex<T>, std::complex<T>, BinaryOp> {
    using ReturnType = stan::math::complex<T>;
  };

  template <typename T, typename BinaryOp>
  struct ScalarBinaryOpTraits<std::complex<T>, stan::math::complex<T>, BinaryOp> {
    using ReturnType = stan::math::complex<T>;
  };


}  // namespace Eigen


namespace Eigen {

namespace internal {
template <typename MatrixType, bool IsComplex>
struct complex_schur_patched_reduce_to_hessenberg;
}

/** \eigenvalues_module \ingroup Eigenvalues_Module
 *
 *
 * \class ComplexSchurPatched
 *
 * \brief Performs a complex Schur decomposition of a real or complex square matrix
 *
 * \tparam MatrixType_ the type of the matrix of which we are
 * computing the Schur decomposition; this is expected to be an
 * instantiation of the Matrix class template.
 *
 * Given a real or complex square matrix A, this class computes the
 * Schur decomposition: \f$ A = U T U^*\f$ where U is a unitary
 * complex matrix, and T is a complex upper triangular matrix.  The
 * diagonal of the matrix T corresponds to the eigenvalues of the
 * matrix A.
 *
 * Call the function compute() to compute the Schur decomposition of
 * a given matrix. Alternatively, you can use the
 * ComplexSchurPatched(const MatrixType&, bool) constructor which computes
 * the Schur decomposition at construction time. Once the
 * decomposition is computed, you can use the matrixU() and matrixT()
 * functions to retrieve the matrices U and V in the decomposition.
 *
 * \note This code is inspired from Jampack
 *
 * \sa class RealSchur, class EigenSolver, class ComplexEigenSolver
 */
template <typename MatrixType_>
class ComplexSchurPatched {
 public:
  typedef MatrixType_ MatrixType;
  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime,
    Options = internal::traits<MatrixType>::Options,
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
  };

  /** \brief Scalar type for matrices of type \p MatrixType_. */
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Eigen::Index Index;  ///< \deprecated since Eigen 3.3

  /** \brief Complex scalar type for \p MatrixType_.
   *
   * This is \c std::complex<Scalar> if #Scalar is real (e.g.,
   * \c float or \c double) and just \c Scalar if #Scalar is
   * complex.
   */
  typedef Scalar ComplexScalar;

  /** \brief Type for the matrices in the Schur decomposition.
   *
   * This is a square matrix with entries of type #ComplexScalar.
   * The size is the same as the size of \p MatrixType_.
   */
  typedef Matrix<ComplexScalar, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime,
                 MaxColsAtCompileTime>
      ComplexMatrixType;

  /** \brief Default constructor.
   *
   * \param [in] size  Positive integer, size of the matrix whose Schur decomposition will be computed.
   *
   * The default constructor is useful in cases in which the user
   * intends to perform decompositions via compute().  The \p size
   * parameter is only used as a hint. It is not an error to give a
   * wrong \p size, but it may impair performance.
   *
   * \sa compute() for an example.
   */
  explicit ComplexSchurPatched(Index size = RowsAtCompileTime == Dynamic ? 1 : RowsAtCompileTime)
      : m_matT(size, size),
        m_matU(size, size),
        m_hess(size),
        m_isInitialized(false),
        m_matUisUptodate(false),
        m_maxIters(-1) {}

  /** \brief Constructor; computes Schur decomposition of given matrix.
   *
   * \param[in]  matrix    Square matrix whose Schur decomposition is to be computed.
   * \param[in]  computeU  If true, both T and U are computed; if false, only T is computed.
   *
   * This constructor calls compute() to compute the Schur decomposition.
   *
   * \sa matrixT() and matrixU() for examples.
   */
  template <typename InputType>
  explicit ComplexSchurPatched(const EigenBase<InputType>& matrix, bool computeU = true)
      : m_matT(matrix.rows(), matrix.cols()),
        m_matU(matrix.rows(), matrix.cols()),
        m_hess(matrix.rows()),
        m_isInitialized(false),
        m_matUisUptodate(false),
        m_maxIters(-1) {
    compute(matrix.derived(), computeU);
  }

  /** \brief Returns the unitary matrix in the Schur decomposition.
   *
   * \returns A const reference to the matrix U.
   *
   * It is assumed that either the constructor
   * ComplexSchurPatched(const MatrixType& matrix, bool computeU) or the
   * member function compute(const MatrixType& matrix, bool computeU)
   * has been called before to compute the Schur decomposition of a
   * matrix, and that \p computeU was set to true (the default
   * value).
   *
   * Example: \include ComplexSchurPatched_matrixU.cpp
   * Output: \verbinclude ComplexSchurPatched_matrixU.out
   */
  const ComplexMatrixType& matrixU() const {
    eigen_assert(m_isInitialized && "ComplexSchurPatched is not initialized.");
    eigen_assert(m_matUisUptodate && "The matrix U has not been computed during the ComplexSchurPatched decomposition.");
    return m_matU;
  }

  /** \brief Returns the triangular matrix in the Schur decomposition.
   *
   * \returns A const reference to the matrix T.
   *
   * It is assumed that either the constructor
   * ComplexSchurPatched(const MatrixType& matrix, bool computeU) or the
   * member function compute(const MatrixType& matrix, bool computeU)
   * has been called before to compute the Schur decomposition of a
   * matrix.
   *
   * Note that this function returns a plain square matrix. If you want to reference
   * only the upper triangular part, use:
   * \code schur.matrixT().triangularView<Upper>() \endcode
   *
   * Example: \include ComplexSchurPatched_matrixT.cpp
   * Output: \verbinclude ComplexSchurPatched_matrixT.out
   */
  const ComplexMatrixType& matrixT() const {
    eigen_assert(m_isInitialized && "ComplexSchurPatched is not initialized.");
    return m_matT;
  }

  /** \brief Computes Schur decomposition of given matrix.
    *
    * \param[in]  matrix  Square matrix whose Schur decomposition is to be computed.
    * \param[in]  computeU  If true, both T and U are computed; if false, only T is computed.

    * \returns    Reference to \c *this
    *
    * The Schur decomposition is computed by first reducing the
    * matrix to Hessenberg form using the class
    * HessenbergDecomposition. The Hessenberg matrix is then reduced
    * to triangular form by performing QR iterations with a single
    * shift. The cost of computing the Schur decomposition depends
    * on the number of iterations; as a rough guide, it may be taken
    * on the number of iterations; as a rough guide, it may be taken
    * to be \f$25n^3\f$ complex flops, or \f$10n^3\f$ complex flops
    * if \a computeU is false.
    *
    * Example: \include ComplexSchurPatched_compute.cpp
    * Output: \verbinclude ComplexSchurPatched_compute.out
    *
    * \sa compute(const MatrixType&, bool, Index)
    */
  template <typename InputType>
  ComplexSchurPatched& compute(const EigenBase<InputType>& matrix, bool computeU = true);

  /** \brief Compute Schur decomposition from a given Hessenberg matrix
   *  \param[in] matrixH Matrix in Hessenberg form H
   *  \param[in] matrixQ orthogonal matrix Q that transform a matrix A to H : A = Q H Q^T
   *  \param computeU Computes the matriX U of the Schur vectors
   * \return Reference to \c *this
   *
   *  This routine assumes that the matrix is already reduced in Hessenberg form matrixH
   *  using either the class HessenbergDecomposition or another mean.
   *  It computes the upper quasi-triangular matrix T of the Schur decomposition of H
   *  When computeU is true, this routine computes the matrix U such that
   *  A = U T U^T =  (QZ) T (QZ)^T = Q H Q^T where A is the initial matrix
   *
   * NOTE Q is referenced if computeU is true; so, if the initial orthogonal matrix
   * is not available, the user should give an identity matrix (Q.setIdentity())
   *
   * \sa compute(const MatrixType&, bool)
   */
  template <typename HessMatrixType, typename OrthMatrixType>
  ComplexSchurPatched& computeFromHessenberg(const HessMatrixType& matrixH, const OrthMatrixType& matrixQ,
                                      bool computeU = true);

  /** \brief Reports whether previous computation was successful.
   *
   * \returns \c Success if computation was successful, \c NoConvergence otherwise.
   */
  ComputationInfo info() const {
    eigen_assert(m_isInitialized && "ComplexSchurPatched is not initialized.");
    return m_info;
  }

  /** \brief Sets the maximum number of iterations allowed.
   *
   * If not specified by the user, the maximum number of iterations is m_maxIterationsPerRow times the size
   * of the matrix.
   */
  ComplexSchurPatched& setMaxIterations(Index maxIters) {
    m_maxIters = maxIters;
    return *this;
  }

  /** \brief Returns the maximum number of iterations. */
  Index getMaxIterations() { return m_maxIters; }

  /** \brief Maximum number of iterations per row.
   *
   * If not otherwise specified, the maximum number of iterations is this number times the size of the
   * matrix. It is currently set to 30.
   */
  static const int m_maxIterationsPerRow = 30;

 protected:
  ComplexMatrixType m_matT, m_matU;
  HessenbergDecomposition<MatrixType> m_hess;
  ComputationInfo m_info;
  bool m_isInitialized;
  bool m_matUisUptodate;
  Index m_maxIters;

 private:
  bool subdiagonalEntryIsNeglegible(Index i);
  ComplexScalar computeShift(Index iu, Index iter);
  void reduceToTriangularForm(bool computeU);
  friend struct internal::complex_schur_patched_reduce_to_hessenberg<MatrixType, NumTraits<Scalar>::IsComplex>;
};

/** If m_matT(i+1,i) is negligible in floating point arithmetic
 * compared to m_matT(i,i) and m_matT(j,j), then set it to zero and
 * return true, else return false. */
template <typename MatrixType>
inline bool ComplexSchurPatched<MatrixType>::subdiagonalEntryIsNeglegible(Index i) {
  RealScalar d = numext::norm1(m_matT.coeff(i, i)) + numext::norm1(m_matT.coeff(i + 1, i + 1));
  RealScalar sd = numext::norm1(m_matT.coeff(i + 1, i));
  if (internal::isMuchSmallerThan(sd, d, NumTraits<RealScalar>::epsilon())) {
    m_matT.coeffRef(i + 1, i) = ComplexScalar(0);
    return true;
  }
  return false;
}

/** Compute the shift in the current QR iteration. */
template <typename MatrixType>
typename ComplexSchurPatched<MatrixType>::ComplexScalar ComplexSchurPatched<MatrixType>::computeShift(Index iu, Index iter) {
  using std::abs;
  if ((iter == 10 || iter == 20) && iu > 1) {
    // exceptional shift, taken from http://www.netlib.org/eispack/comqr.f
    return abs(m_matT.coeff(iu, iu - 1).real()) + abs(m_matT.coeff(iu - 1, iu - 2).real());
  }

  // compute the shift as one of the eigenvalues of t, the 2x2
  // diagonal block on the bottom of the active submatrix
  Matrix<ComplexScalar, 2, 2> t = m_matT.template block<2, 2>(iu - 1, iu - 1);
  RealScalar normt = t.cwiseAbs().sum();
  t /= normt;  // the normalization by sf is to avoid under/overflow

  ComplexScalar b = t.coeff(0, 1) * t.coeff(1, 0);
  ComplexScalar c = t.coeff(0, 0) - t.coeff(1, 1);
  ComplexScalar disc = sqrt(c * c + RealScalar(4) * b);
  ComplexScalar det = t.coeff(0, 0) * t.coeff(1, 1) - b;
  ComplexScalar trace = t.coeff(0, 0) + t.coeff(1, 1);
  ComplexScalar eival1 = (trace + disc) / RealScalar(2);
  ComplexScalar eival2 = (trace - disc) / RealScalar(2);
  RealScalar eival1_norm = numext::norm1(eival1);
  RealScalar eival2_norm = numext::norm1(eival2);
  // A division by zero can only occur if eival1==eival2==0.
  // In this case, det==0, and all we have to do is checking that eival2_norm!=0
  if (eival1_norm > eival2_norm)
    eival2 = det / eival1;
  else if (eival2_norm != RealScalar(0))
    eival1 = det / eival2;

  // choose the eigenvalue closest to the bottom entry of the diagonal
  if (Eigen::numext::norm1(eival1 - t.coeff(1, 1)) < numext::norm1(eival2 - t.coeff(1, 1)))
    return normt * eival1;
  else
    return normt * eival2;
}

template <typename MatrixType>
template <typename InputType>
ComplexSchurPatched<MatrixType>& ComplexSchurPatched<MatrixType>::compute(const EigenBase<InputType>& matrix, bool computeU) {
  m_matUisUptodate = false;
  eigen_assert(matrix.cols() == matrix.rows());

  if (matrix.cols() == 1) {
    m_matT = matrix.derived().template cast<ComplexScalar>();
    if (computeU) m_matU = ComplexMatrixType::Identity(1, 1);
    m_info = Success;
    m_isInitialized = true;
    m_matUisUptodate = computeU;
    return *this;
  }

  internal::complex_schur_patched_reduce_to_hessenberg<MatrixType, NumTraits<Scalar>::IsComplex>::run(*this, matrix.derived(),
                                                                                              computeU);
  computeFromHessenberg(m_matT, m_matU, computeU);
  return *this;
}

template <typename MatrixType>
template <typename HessMatrixType, typename OrthMatrixType>
ComplexSchurPatched<MatrixType>& ComplexSchurPatched<MatrixType>::computeFromHessenberg(const HessMatrixType& matrixH,
                                                                          const OrthMatrixType& matrixQ,
                                                                          bool computeU) {
  m_matT = matrixH;
  if (computeU) m_matU = matrixQ;
  reduceToTriangularForm(computeU);
  return *this;
}
namespace internal {

/* Reduce given matrix to Hessenberg form */
template <typename MatrixType, bool IsComplex>
struct complex_schur_patched_reduce_to_hessenberg {
  // this is the implementation for the case IsComplex = true
  static void run(ComplexSchurPatched<MatrixType>& _this, const MatrixType& matrix, bool computeU) {
    _this.m_hess.compute(matrix);
    _this.m_matT = _this.m_hess.matrixH();
    if (computeU) _this.m_matU = _this.m_hess.matrixQ();
  }
};

template <typename MatrixType>
struct complex_schur_patched_reduce_to_hessenberg<MatrixType, false> {
  static void run(ComplexSchurPatched<MatrixType>& _this, const MatrixType& matrix, bool computeU) {
    typedef typename ComplexSchurPatched<MatrixType>::ComplexScalar ComplexScalar;

    // Note: m_hess is over RealScalar; m_matT and m_matU is over ComplexScalar
    _this.m_hess.compute(matrix);
    _this.m_matT = _this.m_hess.matrixH().template cast<ComplexScalar>();
    if (computeU) {
      // This may cause an allocation which seems to be avoidable
      MatrixType Q = _this.m_hess.matrixQ();
      _this.m_matU = Q.template cast<ComplexScalar>();
    }
  }
};

}  // end namespace internal

// Reduce the Hessenberg matrix m_matT to triangular form by QR iteration.
template <typename MatrixType>
void ComplexSchurPatched<MatrixType>::reduceToTriangularForm(bool computeU) {
  Index maxIters = m_maxIters;
  if (maxIters == -1) maxIters = m_maxIterationsPerRow * m_matT.rows();

  // The matrix m_matT is divided in three parts.
  // Rows 0,...,il-1 are decoupled from the rest because m_matT(il,il-1) is zero.
  // Rows il,...,iu is the part we are working on (the active submatrix).
  // Rows iu+1,...,end are already brought in triangular form.
  Index iu = m_matT.cols() - 1;
  Index il;
  Index iter = 0;       // number of iterations we are working on the (iu,iu) element
  Index totalIter = 0;  // number of iterations for whole matrix

  while (true) {
    // find iu, the bottom row of the active submatrix
    while (iu > 0) {
      if (!subdiagonalEntryIsNeglegible(iu - 1)) break;
      iter = 0;
      --iu;
    }

    // if iu is zero then we are done; the whole matrix is triangularized
    if (iu == 0) break;

    // if we spent too many iterations, we give up
    iter++;
    totalIter++;
    if (totalIter > maxIters) break;

    // find il, the top row of the active submatrix
    il = iu - 1;
    while (il > 0 && !subdiagonalEntryIsNeglegible(il - 1)) {
      --il;
    }

    /* perform the QR step using Givens rotations. The first rotation
       creates a bulge; the (il+2,il) element becomes nonzero. This
       bulge is chased down to the bottom of the active submatrix. */

    ComplexScalar shift = computeShift(iu, iter);
    JacobiRotation<ComplexScalar> rot;
    rot.makeGivens(m_matT.coeff(il, il) - shift, m_matT.coeff(il + 1, il));
    m_matT.rightCols(m_matT.cols() - il).applyOnTheLeft(il, il + 1, rot.adjoint());
    m_matT.topRows((std::min)(il + 2, iu) + 1).applyOnTheRight(il, il + 1, rot);
    if (computeU) m_matU.applyOnTheRight(il, il + 1, rot);

    for (Index i = il + 1; i < iu; i++) {
      rot.makeGivens(m_matT.coeffRef(i, i - 1), m_matT.coeffRef(i + 1, i - 1), &m_matT.coeffRef(i, i - 1));
      m_matT.coeffRef(i + 1, i - 1) = ComplexScalar(0);
      m_matT.rightCols(m_matT.cols() - i).applyOnTheLeft(i, i + 1, rot.adjoint());
      m_matT.topRows((std::min)(i + 2, iu) + 1).applyOnTheRight(i, i + 1, rot);
      if (computeU) m_matU.applyOnTheRight(i, i + 1, rot);
    }
  }

  if (totalIter <= maxIters)
    m_info = Success;
  else
    m_info = NoConvergence;

  m_isInitialized = true;
  m_matUisUptodate = computeU;
}

}  // end namespace Eigen




namespace Eigen {



/** \eigenvalues_module \ingroup Eigenvalues_Module
 *
 *
 * \class ComplexEigenSolverPatched
 *
 * \brief Computes eigenvalues and eigenvectors of general complex matrices
 *
 * \tparam MatrixType_ the type of the matrix of which we are
 * computing the eigendecomposition; this is expected to be an
 * instantiation of the Matrix class template.
 *
 * The eigenvalues and eigenvectors of a matrix \f$ A \f$ are scalars
 * \f$ \lambda \f$ and vectors \f$ v \f$ such that \f$ Av = \lambda v
 * \f$.  If \f$ D \f$ is a diagonal matrix with the eigenvalues on
 * the diagonal, and \f$ V \f$ is a matrix with the eigenvectors as
 * its columns, then \f$ A V = V D \f$. The matrix \f$ V \f$ is
 * almost always invertible, in which case we have \f$ A = V D V^{-1}
 * \f$. This is called the eigendecomposition.
 *
 * The main function in this class is compute(), which computes the
 * eigenvalues and eigenvectors of a given function. The
 * documentation for that function contains an example showing the
 * main features of the class.
 *
 * \sa class EigenSolver, class SelfAdjointEigenSolver
 */
template <typename MatrixType_>
class ComplexEigenSolverPatched {
 public:
  /** \brief Synonym for the template parameter \p MatrixType_. */
  typedef MatrixType_ MatrixType;

  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime,
    Options = internal::traits<MatrixType>::Options,
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
  };

  /** \brief Scalar type for matrices of type #MatrixType. */
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Eigen::Index Index;  ///< \deprecated since Eigen 3.3

  /** \brief Complex scalar type for #MatrixType.
   *
   * This is \c std::complex<Scalar> if #Scalar is real (e.g.,
   * \c float or \c double) and just \c Scalar if #Scalar is
   * complex.
   */
  using ComplexScalar = Scalar;

  /** \brief Type for vector of eigenvalues as returned by eigenvalues().
   *
   * This is a column vector with entries of type #ComplexScalar.
   * The length of the vector is the size of #MatrixType.
   */
  typedef Matrix<ComplexScalar, ColsAtCompileTime, 1, Options & (~RowMajor), MaxColsAtCompileTime, 1> EigenvalueType;

  /** \brief Type for matrix of eigenvectors as returned by eigenvectors().
   *
   * This is a square matrix with entries of type #ComplexScalar.
   * The size is the same as the size of #MatrixType.
   */
  typedef Matrix<ComplexScalar, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime,
                 MaxColsAtCompileTime>
      EigenvectorType;

  /** \brief Default constructor.
   *
   * The default constructor is useful in cases in which the user intends to
   * perform decompositions via compute().
   */
  ComplexEigenSolverPatched()
      : m_eivec(), m_eivalues(), m_schur(), m_isInitialized(false), m_eigenvectorsOk(false), m_matX() {}

  /** \brief Default Constructor with memory preallocation
   *
   * Like the default constructor but with preallocation of the internal data
   * according to the specified problem \a size.
   * \sa ComplexEigenSolverPatched()
   */
  explicit ComplexEigenSolverPatched(Index size)
      : m_eivec(size, size),
        m_eivalues(size),
        m_schur(size),
        m_isInitialized(false),
        m_eigenvectorsOk(false),
        m_matX(size, size) {}

  /** \brief Constructor; computes eigendecomposition of given matrix.
   *
   * \param[in]  matrix  Square matrix whose eigendecomposition is to be computed.
   * \param[in]  computeEigenvectors  If true, both the eigenvectors and the
   *    eigenvalues are computed; if false, only the eigenvalues are
   *    computed.
   *
   * This constructor calls compute() to compute the eigendecomposition.
   */
  template <typename InputType>
  explicit ComplexEigenSolverPatched(const EigenBase<InputType>& matrix, bool computeEigenvectors = true)
      : m_eivec(matrix.rows(), matrix.cols()),
        m_eivalues(matrix.cols()),
        m_schur(matrix.rows()),
        m_isInitialized(false),
        m_eigenvectorsOk(false),
        m_matX(matrix.rows(), matrix.cols()) {
    compute(matrix.derived(), computeEigenvectors);
  }

  /** \brief Returns the eigenvectors of given matrix.
   *
   * \returns  A const reference to the matrix whose columns are the eigenvectors.
   *
   * \pre Either the constructor
   * ComplexEigenSolverPatched(const MatrixType& matrix, bool) or the member
   * function compute(const MatrixType& matrix, bool) has been called before
   * to compute the eigendecomposition of a matrix, and
   * \p computeEigenvectors was set to true (the default).
   *
   * This function returns a matrix whose columns are the eigenvectors. Column
   * \f$ k \f$ is an eigenvector corresponding to eigenvalue number \f$ k
   * \f$ as returned by eigenvalues().  The eigenvectors are normalized to
   * have (Euclidean) norm equal to one. The matrix returned by this
   * function is the matrix \f$ V \f$ in the eigendecomposition \f$ A = V D
   * V^{-1} \f$, if it exists.
   *
   * Example: \include ComplexEigenSolver_eigenvectors.cpp
   * Output: \verbinclude ComplexEigenSolver_eigenvectors.out
   */
  const EigenvectorType& eigenvectors() const {
    eigen_assert(m_isInitialized && "ComplexEigenSolverPatched is not initialized.");
    eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
    return m_eivec;
  }

  /** \brief Returns the eigenvalues of given matrix.
   *
   * \returns A const reference to the column vector containing the eigenvalues.
   *
   * \pre Either the constructor
   * ComplexEigenSolverPatched(const MatrixType& matrix, bool) or the member
   * function compute(const MatrixType& matrix, bool) has been called before
   * to compute the eigendecomposition of a matrix.
   *
   * This function returns a column vector containing the
   * eigenvalues. Eigenvalues are repeated according to their
   * algebraic multiplicity, so there are as many eigenvalues as
   * rows in the matrix. The eigenvalues are not sorted in any particular
   * order.
   *
   * Example: \include ComplexEigenSolver_eigenvalues.cpp
   * Output: \verbinclude ComplexEigenSolver_eigenvalues.out
   */
  const EigenvalueType& eigenvalues() const {
    eigen_assert(m_isInitialized && "ComplexEigenSolverPatched is not initialized.");
    return m_eivalues;
  }

  /** \brief Computes eigendecomposition of given matrix.
   *
   * \param[in]  matrix  Square matrix whose eigendecomposition is to be computed.
   * \param[in]  computeEigenvectors  If true, both the eigenvectors and the
   *    eigenvalues are computed; if false, only the eigenvalues are
   *    computed.
   * \returns    Reference to \c *this
   *
   * This function computes the eigenvalues of the complex matrix \p matrix.
   * The eigenvalues() function can be used to retrieve them.  If
   * \p computeEigenvectors is true, then the eigenvectors are also computed
   * and can be retrieved by calling eigenvectors().
   *
   * The matrix is first reduced to Schur form using the
   * ComplexSchurPatched class. The Schur decomposition is then used to
   * compute the eigenvalues and eigenvectors.
   *
   * The cost of the computation is dominated by the cost of the
   * Schur decomposition, which is \f$ O(n^3) \f$ where \f$ n \f$
   * is the size of the matrix.
   *
   * Example: \include ComplexEigenSolver_compute.cpp
   * Output: \verbinclude ComplexEigenSolver_compute.out
   */
  template <typename InputType>
  ComplexEigenSolverPatched& compute(const EigenBase<InputType>& matrix, bool computeEigenvectors = true);

  /** \brief Reports whether previous computation was successful.
   *
   * \returns \c Success if computation was successful, \c NoConvergence otherwise.
   */
  ComputationInfo info() const {
    eigen_assert(m_isInitialized && "ComplexEigenSolverPatched is not initialized.");
    return m_schur.info();
  }

  /** \brief Sets the maximum number of iterations allowed. */
  ComplexEigenSolverPatched& setMaxIterations(Index maxIters) {
    m_schur.setMaxIterations(maxIters);
    return *this;
  }

  /** \brief Returns the maximum number of iterations. */
  Index getMaxIterations() { return m_schur.getMaxIterations(); }

 protected:
  EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar)

  EigenvectorType m_eivec;
  EigenvalueType m_eivalues;
  ComplexSchurPatched<MatrixType> m_schur;
  bool m_isInitialized;
  bool m_eigenvectorsOk;
  EigenvectorType m_matX;

 private:
  void doComputeEigenvectors(RealScalar matrixnorm);
  void sortEigenvalues(bool computeEigenvectors);
};

template <typename MatrixType>
template <typename InputType>
ComplexEigenSolverPatched<MatrixType>& ComplexEigenSolverPatched<MatrixType>::compute(const EigenBase<InputType>& matrix,
                                                                        bool computeEigenvectors) {
  // this code is inspired from Jampack
  eigen_assert(matrix.cols() == matrix.rows());

  // Do a complex Schur decomposition, A = U T U^*
  // The eigenvalues are on the diagonal of T.
  m_schur.compute(matrix.derived(), computeEigenvectors);

  if (m_schur.info() == Success) {
    m_eivalues = m_schur.matrixT().diagonal();
    if (computeEigenvectors) doComputeEigenvectors(m_schur.matrixT().norm());
    sortEigenvalues(computeEigenvectors);
  }

  m_isInitialized = true;
  m_eigenvectorsOk = computeEigenvectors;
  return *this;
}

template <typename MatrixType>
void ComplexEigenSolverPatched<MatrixType>::doComputeEigenvectors(RealScalar matrixnorm) {
  const Index n = m_eivalues.size();

  matrixnorm = numext::maxi(matrixnorm, (std::numeric_limits<RealScalar>::min)());

  // Compute X such that T = X D X^(-1), where D is the diagonal of T.
  // The matrix X is unit triangular.
  m_matX = EigenvectorType::Zero(n, n);
  for (Index k = n - 1; k >= 0; k--) {
    m_matX.coeffRef(k, k) = ComplexScalar(1.0, 0.0);
    // Compute X(i,k) using the (i,k) entry of the equation X T = D X
    for (Index i = k - 1; i >= 0; i--) {
      m_matX.coeffRef(i, k) = -m_schur.matrixT().coeff(i, k);
      if (k - i - 1 > 0)
        m_matX.coeffRef(i, k) -=
            (m_schur.matrixT().row(i).segment(i + 1, k - i - 1) * m_matX.col(k).segment(i + 1, k - i - 1)).value();
      ComplexScalar z = m_schur.matrixT().coeff(i, i) - m_schur.matrixT().coeff(k, k);
      if (z == ComplexScalar(0)) {
        // If the i-th and k-th eigenvalue are equal, then z equals 0.
        // Use a small value instead, to prevent division by zero.
        numext::real_ref(z) = NumTraits<RealScalar>::epsilon() * matrixnorm;
      }
      m_matX.coeffRef(i, k) = m_matX.coeff(i, k) / z;
    }
  }

  // Compute V as V = U X; now A = U T U^* = U X D X^(-1) U^* = V D V^(-1)
  m_eivec.noalias() = m_schur.matrixU() * m_matX;
  // .. and normalize the eigenvectors
  for (Index k = 0; k < n; k++) {
    m_eivec.col(k).stableNormalize();
  }
}

template <typename MatrixType>
void ComplexEigenSolverPatched<MatrixType>::sortEigenvalues(bool computeEigenvectors) {
  const Index n = m_eivalues.size();
  for (Index i = 0; i < n; i++) {
    Index k;
    m_eivalues.cwiseAbs().tail(n - i).minCoeff(&k);
    if (k != 0) {
      k += i;
      std::swap(m_eivalues[k], m_eivalues[i]);
      if (computeEigenvectors) m_eivec.col(i).swap(m_eivec.col(k));
    }
  }
}

}  // end namespace Eigen


#endif
