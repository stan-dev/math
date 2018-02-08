#ifndef STAN_MATH_REV_MAT_FUNCTOR_INTEGRATOR_DAE_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_INTEGRATOR_DAE_HPP

#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/mat/functor/idas_forward_system.hpp>
#include <stan/math/rev/mat/functor/idas_integrator.hpp>
#include <algorithm>
#include <ostream>
#include <vector>

namespace stan {
namespace math {
/**
 * Return the solutions for the DAE system with residual
 * specified by functor F,
 * given the specified initial state,
 * initial times, times of desired solution, and parameters and
 * data, writing error and warning messages to the specified
 * stream.
 *
 * @tparam DAE type of DAE system
 * @param[in] f functor for the base ordinary differential equation.
 * @param[in] eq_id array for DAE's variable ID(1 for *
 * derivative variables, 0 for algebraic variables).
 * @param[in] yy0 initial state.
 * @param[in] yp0 initial derivative state.
 * @param[in] t0 initial time.
 * @param[in] ts times of the desired solutions, in strictly
 * increasing order, all greater than the initial time.
 * @param[in] theta parameter
 * @param[in] x_r real data
 * @param[in] x_i int data
 * @param[in] rtol relative tolerance passed to IDAS,
 * recommend <10^3.
 * @param[in] atol absolute tolerance passed to IDAS, problem-dependent.
 * @param[in] max_num_steps maximal number of admissable steps
 * between time-points
 * @param[in] msgs message
 * @return a vector of states, each state being a vector of the
 * same size as the state variable, corresponding to a time in ts.
 */
template <typename F, typename Tyy, typename Typ, typename Tpar>
std::vector<std::vector<typename stan::return_type<Tyy, Typ, Tpar>::type> >
integrate_dae(const F& f, const std::vector<int>& eq_id,
              const std::vector<Tyy>& yy0, const std::vector<Typ>& yp0,
              double t0, const std::vector<double>& ts,
              const std::vector<Tpar>& theta, const std::vector<double>& x_r,
              const std::vector<int>& x_i, const double rtol, const double atol,
              const int64_t max_num_steps = idas_integrator::IDAS_MAX_STEPS,
              std::ostream* msgs = nullptr) {
  stan::math::idas_integrator solver(rtol, atol, max_num_steps);
  stan::math::idas_forward_system<F, Tyy, Typ, Tpar> dae{
      f, eq_id, yy0, yp0, theta, x_r, x_i, msgs};

  return solver.integrate(dae, t0, ts);
}
}  // namespace math
}  // namespace stan

#endif
