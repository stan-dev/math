#ifndef STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LOGIT_LPMF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LOGIT_LPMF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#include <string>

namespace stan {
  namespace math {

    /**
     * Returns the log PMF of the logit-parametrized Bernoulli distribution. If
     * containers are supplied, returns the log sum of the probabilities.
     *
     * @tparam T_n type of integer parameter
     * @tparam T_prob type of chance of success parameter
     * @param n integer parameter
     * @param theta logit-transformed chance of success parameter
     * @return log probability or log sum of probabilities
     * @throw std::domain_error if theta is infinite.
     * @throw std::invalid_argument if container sizes mismatch.
     */
    template <bool propto, typename T_n, typename T_prob> // propto: false means we care about the normalisation constant.
    typename return_type<T_prob>::type
    bernoulli_logit_lpmf(const T_n& n, const T_prob& theta) {
      static const std::string function = "bernoulli_logit_lpmf";
      typedef typename stan::partials_return_type<T_n, T_prob>::type
        T_partials_return; // One level down: var-> double

      using stan::is_constant_struct;
      using std::exp;

      if (!(stan::length(n) && stan::length(theta)))
        return 0.0;

      T_partials_return logp(0.0); // This initialises log probability to 0.0.
//some error checking next
      check_bounded(function, "n", n, 0, 1);  // Right, so n is supposed to hold 0s and 1s (it is our usual y)
      check_not_nan(function, "Logit transformed probability parameter", theta);
      check_consistent_sizes(function,
                             "Random variable", n,
                             "Probability parameter", theta);

      if (!include_summand<propto, T_prob>::value)
        return 0.0;

      scalar_seq_view<T_n> n_vec(n); //scalar seq view takes scalar or sequence and makes it look like a sequence
      scalar_seq_view<T_prob> theta_vec(theta);
      size_t N = max_size(n, theta);
      operands_and_partials<T_prob> ops_partials(theta); // operands_and_partials stores the partial derivatives and creates the autodiff node, the var to return

      for (size_t n = 0; n < N; n++) {
        const int n_int = value_of(n_vec[n]); // value_of gives value of an autodiff variable
        const T_partials_return theta_dbl = value_of(theta_vec[n]);

        const int sign = 2 * n_int - 1; // maps 0,1 to -1,1
        const T_partials_return ntheta = sign * theta_dbl;
        const T_partials_return exp_m_ntheta = exp(-ntheta);

        // Handle extreme values gracefully using Taylor approximations.
        static const double cutoff = 20.0;
        if (ntheta > cutoff)
          logp -= exp_m_ntheta;
        else if (ntheta < -cutoff)
          logp += ntheta;
        else
          logp -= log1p(exp_m_ntheta);

        if (!is_constant_struct<T_prob>::value) {  //is_constant_struct figures out if we just passed in a double (or int), rather than a var
          static const double cutoff = 20.0; // Presumably, we can delete this line.
          if (ntheta > cutoff)
            ops_partials.edge1_.partials_[n] -= exp_m_ntheta; // edge1 has the partial derivatives for the first argument (operand); it always returns a list
          else if (ntheta < -cutoff)
            ops_partials.edge1_.partials_[n] += sign;
          else
            ops_partials.edge1_.partials_[n] += sign * exp_m_ntheta
              / (exp_m_ntheta + 1);
        }
      }
      return ops_partials.build(logp); // ops_partials.build returns a fake autodiff variable (which contains the analytic derivatives) if we start from autodiff 
    }

    template <typename T_n,
              typename T_prob>
    inline
    typename return_type<T_prob>::type
    bernoulli_logit_lpmf(const T_n& n,
                        const T_prob& theta) {
      return bernoulli_logit_lpmf<false>(n, theta);
    }

  }
}
#endif
