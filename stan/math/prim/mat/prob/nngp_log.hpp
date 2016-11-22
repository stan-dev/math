#ifndef NNGP_LOG_HPP
#define NNGP_LOG_HPP

#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/check_greater.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/fun/to_row_vector.hpp>
#include <stan/math/prim/mat/fun/mdivide_left_tri_low.hpp>
#include <stan/math/prim/mat/fun/mdivide_right_tri_low.hpp>
#include <stan/math/prim/mat/fun/square.hpp>

namespace stan {
    namespace math {
        
        /**
         * The log of the Nearest Neighbors Gaussian Process(NNGP)
         * with an exponential kernal for the given nearest 
         * neighbors distance matrics, index of the nearest 
         * neighbors and parameters including partial sill, 
         * nuggets, and range.
         *
         *
         * @param sigma         square root of Partial Sill
         * @param tau           square root of Nuggets
         * @param phi           inverse of Range
         * @param beta          vector of regression coefficients
         * @return              The log of the NNGP density.
         * @tparam T_y, T_x     Type of scalar.
         * @tparam T_beta, T_phi, T_sigma
         *                      Type of parametrs.
         * @tparam T_dist, T_distM
         *                      Type of distance.
         * @tparam T_ind        Type of index of location
         */

        template <typename T_y, typename T_x, typename T_beta,
        typename T_sigma, typename T_tau, typename T_phi,
        typename T_dist, typename T_distM, typename T_ind>
        typename return_type<T_y, T_beta, T_sigma, T_tau, T_phi>::type
        
        nngp_log(const T_y& Y, const T_x& X, const T_beta& beta,
                 const T_sigma& sigma, const T_tau& tau,
                 const T_phi& phi, const T_dist& neardist,
                 const T_distM& neardistM, const T_ind& nearind){
            
            
            static const char* function("nearest_neighbor_GP_log");
            
            typedef typename return_type<T_y, T_beta, T_sigma,
            T_tau, T_phi>::type lp_type;
            lp_type lp(0.0);
            
            int M = neardist.cols();    // no. of neighbors
            int N = X.rows();           // no. of observations
            int P = X.cols();           // no. of predictors
            
            using Eigen::Dynamic;
            using Eigen::Matrix;
            using Eigen::VectorXd;
            
            /* check dimension matching: */
            check_size_match(function,"No. of columes in neardist", M,
                        "No. of columes in nearind", nearind.cols());
            
            /* check positivity: */
            check_greater(function, "phi", phi, 0);
            
            /* check finite and NA*/
            for (size_t i = 0; i < P; i++) {
                check_finite(function,
                             "regression coefficients", beta[i]);
                check_not_nan(function,
                              "regression coefficients", beta[i]);
            }
            check_finite(function, "inverse of Range", phi);
            check_not_nan(function, "inverse of Range", phi);
            check_finite(function, "square root of partial sill",
                         sigma);
            check_not_nan(function, "square root of partial sill",
                          sigma);
            check_finite(function, "square root of nuggets", tau);
            check_not_nan(function, "square root of nuggets", tau);
            
            VectorXd V(N);
            V(0) = (square(tau) + square(sigma))/ square(sigma) ;
            
            VectorXd Uw(N);       // Uw(N) = temp_w(N) = Y - X * beta
            VectorXd temp_w(N);

            Uw = subtract(Y, multiply(X, beta));
            temp_w = Uw;
            
            // get exp(-phi * dist):
            for (int i = 0; i < (N - 1); i++){
                const int dim = (i < M) ? (i + 1) : M;
                Matrix<double, Dynamic, Dynamic> temp_distM(dim, dim);
                int h = 0;
                for (int j = 0; j < dim; j++){
                    for (int k = j; k < dim; k++){
                        temp_distM(j, k) =
                        std::exp(-phi * neardistM(i, j*dim + k));
                        
                        temp_distM(k, j) = temp_distM(j, k);
                    }
                }
                for(int j = 0; j < dim; j++){
                    temp_distM(j, j) =
                    (square(tau) + square(sigma))/ square(sigma);
                }
                
                check_positive(function, "Covariance matrix rows",
                               temp_distM.rows());
                Matrix<double, Dynamic, Dynamic>
                choldistM(cholesky_decompose(temp_distM));

                VectorXd u(dim);
                VectorXd v(dim);
                Matrix<double, 1, Dynamic> v2(1, dim);
                for (int j = 0; j < dim; j++){
                    u(j) = std::exp(-phi * neardist(i, j));
                }

                // v = L^-1 * u

                v = mdivide_left_tri_low(choldistM, u);
                V((i+1)) = (square(tau) + square(sigma))/square(sigma)
                            - stan::math::square(v).sum();
                v2 = mdivide_right_tri_low(to_row_vector(v),
                                           choldistM);
                
                for (int j = 0; j < dim; j++){
                    Uw(i + 1) -= v2(j) * temp_w(nearind(i, j) - 1);
                }
            }
            
            for (int i = 0; i < N; i++){
                lp += -  log(V(i)) -
                1.0 / square(sigma) * (Uw(i) * Uw(i) / V(i));
            }
            
            lp = 0.5 * (lp - N * log(square(sigma)));
            
            return lp;
        }
        
    }
}
#endif
















