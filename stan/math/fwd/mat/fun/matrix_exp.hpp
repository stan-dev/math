#ifndef STAN_MATH_FWD_MAT_FUN_MATRIX_EXP_HPP
#define STAN_MATH_FWD_MAT_FUN_MATRIX_EXP_HPP

#include <unsupported/Eigen/MatrixFunctions>
#include <stan/math/prim/mat/fun/MatrixExponential.h>
#include <stan/math/prim/mat/fun/matrix_exp_spec.hpp>


namespace stan {
    namespace math {
        
        using Eigen::Matrix;
        using Eigen::Dynamic;
          
        template <typename MatrixType, typename T>
        struct matrix_exp_computeUV<MatrixType, fvar<T> >
        {
          static void run(const MatrixType& arg, MatrixType& U, MatrixType& V, int& squarings)
          {
            using std::frexp;
            using std::pow;
            const fvar<T> l1norm = arg.cwiseAbs().colwise().sum().maxCoeff();
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
          * Computes the derivative of matrix raised to a power
          * with respect to a parameter, by recursively applying 
          * the chain rule. 
          *
          * @param A input matrix
          * @param k the exponent of the matrix
          * @return matrix of derivatives
          */
        
        template<typename T>
        Matrix<fvar<T>, Dynamic, Dynamic> 
        dpow(Matrix<fvar<T>, Dynamic, Dynamic>, k) {  
        
        	
        	
        
        
        /**
          * Computes the derivative of the matrix exponential, 
          * obtained through a Pade approximation.
          *
          * @param A input matrix
          * @param degree of Pade approximation (assume p = q).
          */
        
        template<typename T>
        void
        dexpm(const Matrix<fvar<T>, Dynamic, Dynamic>& A, int p) {
        
        	Matrix<T, Dynamic, Dynamic> dA(A.rows(), A.cols());
        	Matrix<T, Dynamic, Dynamic> dN(A.rows(), A.cols())
        		= Matrix<T, Dynamic, Dynamic>::Zero(A.rows(), A.cols());
        	Matrix<T, Dynamic, Dynamic> dD(A.rows(), A.cols())
        		= Matrix<T, Dynamic, Dynamic>::Zero(A.rows(), A.cols());
        	
        	for (int i = 0; i < A.rows(); i++) {
        		for (int j = 0; i < A.cols(); i++) {
        			dA(i,j) = A(i,j).d_;
        		}
        	}
        	
        	for (int i = 0; i <= p; i++) {
        		dN += factorial(2*p - i) * factorial(p) * dpow(A, i)
        				/ (factorial(2*p) * factorial(i) * factorial(p-k));
        				
        		dD += factorial(2*p - i) * factorial(q) * 
        		
        	
        
        	
        	
        
        	          
          
        
        
       /**
         * Returns the exponential of a matrix. 
         *
         * @param A input matrix
         * @return B matrix exponential of A
         */
      
       template<typename T>
       inline
       Eigen::Matrix<fvar<T>, Eigen::Dynamic, Eigen::Dynamic>
       expm(const Eigen::Matrix<fvar<T>, Eigen::Dynamic, Eigen::Dynamic>& A) {
    
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
