#ifndef STAN_MATH_PRIM_ARR_FUN_DECOUPLE_ODE_STATES_HPP
#define STAN_MATH_PRIM_ARR_FUN_DECOUPLE_ODE_STATES_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <vector>

namespace stan {

  namespace math {

    /**
     * Takes sensitivity output from integrators and returns results
     * in precomputed_gradients format.
     *
     * Solution input vector size depends on requested sensitivities,
     * which can be enabled for initials and parameters. For each
     * sensitivity N states are computed. The input vector is expected
     * to be ordered by states, i.e. first the N states, then
     * optionally the N sensitivities for the initials (first the N
     * states for the first initial and so on), finally the
     * sensitivities for the M parameters follow optionally.
     *
     * @tparam T1_initial type of scalars for initial values.
     * @tparam T2_param type of scalars for parameters.
     * @param[in] y output from integrator in column-major format
     * as a coupled system output
     * @param[in] y0 initial state.
     * @param[in] theta parameters
     * @return a vector of states for each entry in y in Stan var
     * format
     */
    template <typename T_initial, typename T_param>
    inline
    std::vector<std::vector<typename stan::return_type<T_initial,
                                                       T_param>::type> >
    decouple_ode_states(const std::vector<std::vector<double> >& y,
                        const std::vector<T_initial>& y0,
                        const std::vector<T_param>& theta);

    /**
     * The decouple ODE states operation for the case of no
     * sensitivities is equal to the identity operation.
     *
     * @param[in] y output from integrator 
     * @param[in] y0 initial state.
     * @param[in] theta parameters
     * @return y
     */
    template <>
    inline
    std::vector<std::vector<double> >
    decouple_ode_states(const std::vector<std::vector<double> >& y,
                        const std::vector<double>& y0,
                        const std::vector<double>& theta) {
      return y;
    }

  }
}
#endif
