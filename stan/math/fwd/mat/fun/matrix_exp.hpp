#ifndef STAN_MATH_FWD_MAT_FUN_MATRIX_EXP_HPP
#define STAN_MATH_FWD_MAT_FUN_MATRIX_EXP_HPP

#include <unsupported/Eigen/MatrixFunctions>
#include <stan/math/prim/mat/fun/MatrixExponential.h>
#include <stan/math/prim/mat/fun/matrix_exp_spec.hpp>


namespace stan {
    namespace math {
        
        template <typename MatrixType, typename T>
        struct matrix_exp_computeUV<MatrixType, fvar<T> >
        {
          static void run(const MatrixType& arg, MatrixType& U, MatrixType& V, int& squarings)
          {
            using std::frexp;
            using std::pow;
            const stan::math::fvar<T> l1norm = arg.cwiseAbs().colwise().sum().maxCoeff();
            squarings = 0;
            if (l1norm < 1.495585217958292e-002) {
              matrix_exp_pade3(arg, U, V);
            } else if (l1norm < 2.539398330063230e-001) {
              matrix_exp_pade5(arg, U, V);
            } else if (l1norm < 9.504178996162932e-001) {
              matrix_exp_pade7(arg, U, V);
            } else if (l1norm < 2.097847961257068e+000) {
              matrix_exp_pade9(arg, U, V);
            } else {
              const double maxnorm = 5.371920351148152;
              frexp(l1norm / maxnorm, &squarings);
              if (squarings < 0) squarings = 0;
              MatrixType A = arg.unaryExpr(MatrixExponentialScalingOp<stan::math::fvar<T> >(squarings));
              matrix_exp_pade13(A, U, V);
            }
          }
        };
        
        
       /**
         * Returns the exponential of a matrix. 
         *
         * @param A input matrix
         * @return B matrix exponential of A
         */
      
       template<typename T>
       inline
       Eigen::Matrix<fvar<T>, Eigen::Dynamic, Eigen::Dynamic>
       matrix_exp(const Eigen::Matrix<fvar<T>, Eigen::Dynamic, Eigen::Dynamic>& A) {
    
        check_nonzero_size("matrix_exp", "input matrix", A);
        check_square("matrix_exp", "input matrix", A);
           
        Matrix<fvar<T>, Dynamic, Dynamic> B;
           
        if (A.cols() == 2) matrix_exp_compute_2x2(A, B);
        else if (is_symmetric(A)) matrix_exp_compute_sym(A, B);
        else matrix_exp_compute(A, B);
    
        return B; 
       }

  }
} // end namespace stan::math

 
#endif
