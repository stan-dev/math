#ifndef STAN_MATH_PRIM_MAT_PROB_MULTI_NORMAL_RNG_HPP
#define STAN_MATH_PRIM_MAT_PROB_MULTI_NORMAL_RNG_HPP

#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/mat/err/check_ldlt_factor.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/mat/fun/trace_inv_quad_form_ldlt.hpp>
#include <stan/math/prim/mat/fun/log_determinant_ldlt.hpp>

#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>

namespace stan {
  namespace math {

    /**
     * Return a pseudo-random vector with a multi-variate normal
     * distribution given the specified location parameter and
     * covariance matrix and pseudo-random number generator.
     *
     * If calculating more than one multivariate-normal random draw
     * then it is more efficient to calculate the Cholesky factor of
     * the covariance matrix and use the function
     * <code>stan::math::multi_normal_cholesky_rng</code>.
     *
     * @tparam RNG Type of pseudo-random number generator.
     * @param mu Location parameter.
     * @param S Covariance parameter.
     * @param rng Pseudo-random number generator.
     */
    template <class RNG>
    inline Eigen::VectorXd
    multi_normal_rng(const Eigen::VectorXd& mu, const Eigen::MatrixXd& S,
                     RNG& rng) {
      using boost::variate_generator;
      using boost::normal_distribution;

      static const char* function("multi_normal_rng");

      check_positive(function, "Covariance matrix rows", S.rows());
      check_symmetric(function, "Covariance matrix", S);
      check_finite(function, "Location parameter", mu);

      Eigen::LLT<Eigen::MatrixXd> llt_of_S = S.llt();
      if (llt_of_S.info() != Eigen::Success) {
        throw std::domain_error("multi_normal_rng:"
                                " Matrix argument not solvable;"
                                " probably not positive definite.");
      }

      variate_generator<RNG&, normal_distribution<> >
        std_normal_rng(rng, normal_distribution<>(0, 1));

      Eigen::VectorXd z(S.cols());
      for (int i = 0; i < S.cols(); i++)
        z(i) = std_normal_rng();

      return mu + llt_of_S.matrixL() * z;
    }

  }
}
#endif
