#ifndef STAN_MATH_REV_ARR_FUNCTOR_DECOUPLE_STATES_HPP
#define STAN_MATH_REV_ARR_FUNCTOR_DECOUPLE_STATES_HPP

#include <stan/math/rev/core.hpp>

#include <boost/type_traits/is_same.hpp>
#include <vector>

namespace stan {

  namespace math {

    /**
     * Takes sensitivity output from integrators and returns results
     * in precomputed_gradients format.
     *
     * Solution input vector size depends on requested sensitivities,
     * which can be enabled for initials and parameters. Per
     * sensitivity N states are computed. The input vector is expected
     * to be ordered by states, i.e. first the N states, then
     * optionally the N sensitivities for the initials (first the N
     * states for the first initial and so on), finally the
     * sensitivities for the M parameters follow optionally.
     *
     */
    template <typename T1, typename T2>
    inline
    std::vector<std::vector<typename stan::return_type<T1, T2>::type> >
    decouple_states(const std::vector<std::vector<double> >& y,
                    const std::vector<T1>& y0,
                    const std::vector<T2>& theta) {
      using std::vector;
      using stan::math::var;
      using stan::math::precomputed_gradients;

      vector<typename stan::return_type<T1, T2>::type> vars;
      typedef boost::is_same<T1, stan::math::var> initial_var;
      typedef boost::is_same<T2, stan::math::var> theta_var;

      const size_t N = y0.size();
      const size_t M = theta.size();
      const size_t S =
        (initial_var::value ? N : 0) +
        (theta_var::value ? M : 0);

      if (initial_var::value)
        vars.insert(vars.end(), y0.begin(), y0.end());
      if (theta_var::value)
        vars.insert(vars.end(), theta.begin(), theta.end());


      vector<var> temp_vars(N);
      vector<double> temp_gradients(S);
      vector<vector<var> > y_return(y.size());

      for (size_t i = 0; i < y.size(); i++) {
        // iterate over number of equations
        for (size_t j = 0; j < N; j++) {
          // iterate over parameters for each equation
          for (size_t k = 0; k < S; k++)
            temp_gradients[k] = y[i][N + N * k + j];

          temp_vars[j] = precomputed_gradients(y[i][j],
                                               vars, temp_gradients);
        }
        y_return[i] = temp_vars;
      }
      return y_return;
    }

    // special case if both (initials and parameters) are known, then
    // this function just returns its input.
    template <>
    inline
    std::vector<std::vector<double> >
    decouple_states(const std::vector<std::vector<double> >& y,
                    const std::vector<double>& y0,
                    const std::vector<double>& theta) {
      return y;
    }
    
  }  // ns math
}  // ns stan

#endif
