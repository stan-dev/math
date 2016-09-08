#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_2X2_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_2X2_HPP

#include <unsupported/Eigen/MatrixFunctions>
#include <stan/math/prim/scal/err/check_size_match.hpp>
//#include <stan/math/prim/mat/err/check_2x2.hpp>

namespace stan {
    namespace math {   
        
        using Eigen::Matrix;
        using Eigen::Dynamic;
        
        /**
         * Return the matrix exponential of a 2x2 matrix. Reference for 
         * algorithm: http://mathworld.wolfram.com/MatrixExponential.html
         *
         * @param[in] A A 2x2 matrix
         * @return Matrix exponential of A.
         */
        template <typename T>
        Matrix<T, Dynamic, Dynamic>
        matrix_exp_compute_2x2(const Matrix<T, Dynamic, Dynamic>& A) {
            
            Matrix<T, Dynamic, Dynamic> B;
                       
            T a = A(0,0), b = A(0,1), c = A(1,0), d = A(1,1), delta;
            delta = sqrt((a - d)*(a - d) + 4*b*c);
            
            B.resize(2,2);
            B(0,0) = exp(.5*(a+d)) * (delta*cosh(.5*delta) + (a-d)*sinh(.5*delta));
            B(0,1) = 2*b*exp(.5*(a+d)) * sinh(.5*delta);
            B(1,0) = 2*c*exp(.5*(a+d)) * sinh(.5*delta);
            B(1,1) = exp(.5*(a+d)) * (delta*cosh(.5*delta) + (d-a)*sinh(.5*delta));
            
            B = 1/delta * B;
            
            return B;
        }
        
    }
}
#endif