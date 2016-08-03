#ifndef STAN_MATH_REV_MAT_FUN_EIGEN_NUMTRAITS_HPP
#define STAN_MATH_REV_MAT_FUN_EIGEN_NUMTRAITS_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/std_numeric_limits.hpp>
#include <limits>

namespace Eigen {

  /**
   * Numerical traits template override for Eigen for automatic
   * gradient variables.
   *
   * Documentation here: 
   *   http://eigen.tuxfamily.org/dox/structEigen_1_1NumTraits.html
   */
  template<> struct NumTraits<stan::math::var>
    : GenericNumTraits<stan::math::var> {
    static inline stan::math::var dummy_precision() {
      return NumTraits<double>::dummy_precision();
    }

    enum {
      RequireInitialization = 0,
      // ReadCost = twice the cost to copy a double
      ReadCost = 2 * NumTraits<double>::ReadCost,
      // AddCost = single addition going forward
      AddCost = NumTraits<double>::AddCost,
      // MulCost = single multiplication going forward
      MulCost = NumTraits<double>::MulCost
    };
  };

  namespace internal {
    /**
     * Implemented this for printing to stream.
     */
    template<>
    struct significant_decimals_default_impl<stan::math::var, false> {
      static inline int run() {
        using std::ceil;
        using std::log;
        return cast<double, int>(ceil(-log(std::numeric_limits<double>
                                           ::epsilon())
                                      / log(10.0)));
      }
    };

    /**
     * Scalar product traits override for Eigen for automatic
     * gradient variables.
     */
    template <>
    struct scalar_product_traits<stan::math::var, double> {
      typedef stan::math::var ReturnType;
    };

    /**
     * Scalar product traits override for Eigen for automatic
     * gradient variables.
     */
    template <>
    struct scalar_product_traits<double, stan::math::var> {
      typedef stan::math::var ReturnType;
    };

    /**
     * Override matrix-vector and matrix-matrix products to use more efficient implementation.
     */
    template<typename Index, bool ConjugateLhs, bool ConjugateRhs>
    struct general_matrix_vector_product<Index, stan::math::var, ColMajor,
                                         ConjugateLhs, stan::math::var,
                                         ConjugateRhs> {
      typedef stan::math::var LhsScalar;
      typedef stan::math::var RhsScalar;
      typedef typename scalar_product_traits<LhsScalar, RhsScalar>::ReturnType
      ResScalar;
      enum { LhsStorageOrder = ColMajor };

      EIGEN_DONT_INLINE static void run(
                                        Index rows, Index cols,
                                        const LhsScalar* lhs, Index lhsStride,
                                        const RhsScalar* rhs, Index rhsIncr,
                                        ResScalar* res, Index resIncr,
                                        const ResScalar &alpha) {
        for (Index i = 0; i < rows; i++) {
          res[i*resIncr]
            += stan::math::var
            (new stan::math::gevv_vvv_vari
             (&alpha,
              (static_cast<int>(LhsStorageOrder) == static_cast<int>(ColMajor))
              ?(&lhs[i]):(&lhs[i*lhsStride]),
              (static_cast<int>(LhsStorageOrder) == static_cast<int>(ColMajor))
              ?(lhsStride):(1),
              rhs, rhsIncr, cols));
        }
      }
    };
    template<typename Index, bool ConjugateLhs, bool ConjugateRhs>
    struct general_matrix_vector_product<Index, stan::math::var,
                                         RowMajor, ConjugateLhs,
                                         stan::math::var, ConjugateRhs> {
      typedef stan::math::var LhsScalar;
      typedef stan::math::var RhsScalar;
      typedef typename scalar_product_traits<LhsScalar, RhsScalar>::ReturnType
      ResScalar;
      enum { LhsStorageOrder = RowMajor };

      EIGEN_DONT_INLINE static void
      run(Index rows, Index cols,
          const LhsScalar* lhs, Index lhsStride,
          const RhsScalar* rhs, Index rhsIncr,
          ResScalar* res, Index resIncr, const RhsScalar &alpha) {
        for (Index i = 0; i < rows; i++) {
          res[i*resIncr]
            += stan::math::var
            (new stan::math::gevv_vvv_vari
             (&alpha,
              (static_cast<int>(LhsStorageOrder) == static_cast<int>(ColMajor))
              ? (&lhs[i]) : (&lhs[i*lhsStride]),
              (static_cast<int>(LhsStorageOrder) == static_cast<int>(ColMajor))
              ? (lhsStride) : (1),
              rhs, rhsIncr, cols));
        }
      }
    };
    template<typename Index, int LhsStorageOrder, bool ConjugateLhs,
             int RhsStorageOrder, bool ConjugateRhs>
    struct general_matrix_matrix_product<Index, stan::math::var,
                                         LhsStorageOrder, ConjugateLhs,
                                         stan::math::var, RhsStorageOrder,
                                         ConjugateRhs, ColMajor> {
      typedef stan::math::var LhsScalar;
      typedef stan::math::var RhsScalar;
      typedef typename scalar_product_traits<LhsScalar, RhsScalar>::ReturnType
      ResScalar;
      static void run(Index rows, Index cols, Index depth,
                      const LhsScalar* lhs, Index lhsStride,
                      const RhsScalar* rhs, Index rhsStride,
                      ResScalar* res, Index resStride,
                      const ResScalar &alpha,
                      level3_blocking<LhsScalar, RhsScalar>& /* blocking */,
                      GemmParallelInfo<Index>* /* info = 0 */) {
        for (Index i = 0; i < cols; i++) {
          general_matrix_vector_product<Index, LhsScalar, LhsStorageOrder,
                                        ConjugateLhs, RhsScalar, ConjugateRhs>
            ::run(rows, depth, lhs, lhsStride,
                  &rhs[(static_cast<int>(RhsStorageOrder)
                         == static_cast<int>(ColMajor))
                        ? (i*rhsStride) :(i) ],
                  (static_cast<int>(RhsStorageOrder)
                   == static_cast<int>(ColMajor)) ? (1) : (rhsStride),
                  &res[i*resStride], 1, alpha);
        }
      }
    };
  }
}

#endif
