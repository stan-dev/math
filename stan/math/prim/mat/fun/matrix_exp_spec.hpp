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
            
            Matrix<T, Dynamic, Dynamic> landa = diag_matrix(eigenvalues_sym(A)),
                                        V = eigenvectors_sym(A),
                                        V_inv = inverse(V);
            
            B = V * exp(landa); 
            B *= V_inv;
            
        }
        
        /**
         * Return the matrix expontential of a nilpotent matrix.
         *
         * @param A A nilpotent matrix.
         * @return Matrix exponential of A. 
         */
        
        /*template <typename T>
        void matrix_exp_compute_nil(const Matrix<T, Dynamic, Dynamic>& A,
                                    Matrix<T, Dynamic, Dynamic>& B) {
            
            
        }*/
            
            
        
    }
}
#endif