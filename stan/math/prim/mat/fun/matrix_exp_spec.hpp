#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_SPEC_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_SPEC_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
    namespace math {
        
        using Eigen::Matrix;
        using Eigen::Dynamic;
        
        /**
         * Return the matrix exponential for a 2x2 matrix.
         * Source for algorithm: http://mathworld.wolfram.com/MatrixExponential.html
         *
         * @param A A 2x2 matrix
         * @return Matrix exponential of A.
         */
        template <typename T>
        void
        matrix_exp_compute_2x2(const Matrix<T, Dynamic, Dynamic>& A,
                               Matrix<T, Dynamic, Dynamic>& B) {
            
            assert((A.cols() == 2) && (A.rows() == 2));
            
            T a = A(0,0), b = A(0,1), c = A(1,0), d = A(1,1), delta;
            delta = sqrt((a - d)*(a - d) + 4*b*c);
            
            B.resize(2,2);
            B(0,0) = exp(.5*(a+d)) * (delta*cosh(.5*delta) + (a-d)*sinh(.5*delta));
            B(0,1) = 2*b*exp(.5*(a+d)) * sinh(.5*delta);
            B(1,0) = 2*c*exp(.5*(a+d)) * sinh(.5*delta);
            B(1,1) = exp(.5*(a+d)) * (delta*cosh(.5*delta) + (d-a)*sinh(.5*delta));
            
            B = 1/delta * B;
        }
        
        /** 
         * Return the matrix exponential of a symmetric matrix.
         *
         * @param A A symmetric matrix.
         * @return Matrix exponential of A.
         */
        template <typename T>
        void matrix_exp_compute_sym(const Matrix<T, Dynamic, Dynamic>& A,
                                    Matrix<T, Dynamic, Dynamic>& B) {
            
            Matrix<T, Dynamic, Dynamic> landa = diag_matrix(exp(eigenvalues_sym(A))),
                                        V = eigenvectors_sym(A),
                                        V_inv = inverse(V);
            
            B = V * landa;
            B *= V_inv;
        }
        
        
        /**
         * Returns a boolean for whether a matrix is nilpotent 
         * or not. 
         *
         * @input: A a matrix
         * @return: true if the matrix is nilpotent, else false. 
         */ 
         template <typename T> 
         bool is_nilpotent(const Matrix<T, Dynamic, Dynamic> A) {
         	
         	int i, j;
         
         	if (A.rows() == 1) return false;
         	for (i = 0; i < A.rows(); i++) {
         		for (j = 0; j < i; j++) if (A(i,j) != 0) return false; 
         		}
         	return true;
         }
         
        /**
         * Compute the factorial of an integer.
         *
         * @param[in] a integer
         * @return the factorial of a.
         */
        int factorial(int a) {
        	
        	if (a == 0 || a == 1) return 1;
        	else return a * factorial(a - 1);
        }
        
        /**
         * Compute the power of a matrix. 
         *
         * @param[in] A the input matrix.
         * @param[in] a the exponent of the matrix.
         * @return a matrix that is the power of the input matrix.
         */
        template <typename T>
        Eigen::Matrix<T, Dynamic, Dynamic> 
        pow(Eigen::Matrix<T, Dynamic, Dynamic> A, int a) {
                
        	Matrix<T, Dynamic, Dynamic> B 
        				= Matrix<T, Dynamic, Dynamic>::Identity(A.cols(), A.rows());
        	
        	if(a != 0) for (int i = 0; i < a; i++) B *= A;
        	
        	return B;
        }
        	
        /**
         * Return the matrix expontential of a nilpotent matrix.
         *
         * @param A A nilpotent matrix.
         * @return Matrix exponential of A. 
         */
        template <typename T>
        void matrix_exp_compute_nil(const Matrix<T, Dynamic, Dynamic>& A,
                                    Matrix<T, Dynamic, Dynamic>& B) {
                                    
        	int n = A.cols();
        	B = Matrix<T, Dynamic, Dynamic>::Identity(A.cols(), A.rows());
        	
        	for (int i = 1; i < n;  i++) B += (1 / factorial(i)) * pow(A, i);         	
        }                    
    
    } // end namespace math
} // end namespace stan

#endif