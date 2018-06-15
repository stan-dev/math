#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_FORWARD_PDE_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_FORWARD_PDE_HPP

#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the solutions for the specified PDE problem at
 * specified location.
 *
 * This function is templated to allow various PDE library
 * interfaces and corresponding input deck.
 *
 * @tparam F_pde type of PDE system interface.
 * @param[in] f functor for the base ordinary differential equation.
 * @param[in] y0 initial state.
 * @param[in] t0 initial time.
 * @param[in] ts times of the desired solutions, in strictly
 * increasing order, all greater than the initial time.
 * @param[in] theta parameter vector for the ODE.
 * @param[in] x continuous data vector for the ODE.
 * @param[in] x_int integer data vector for the ODE.
 * @param[out] msgs the print stream for warning messages.
 * @param[in] relative_tolerance relative tolerance parameter
 *   for Boost's ode solver. Defaults to 1e-6.
 * @param[in] absolute_tolerance absolute tolerance parameter
 *   for Boost's ode solver. Defaults to 1e-6.
 * @param[in] max_num_steps maximum number of steps to take within
 *   the Boost ode solver.
 * @return a vector of states, each state being a vector of the
 * same size as the state variable, corresponding to a time in ts.
 */
  template<typename F_pde>
  inline std::vector<double> forward_pde(const F_pde& pde,
                                         const std::vector<double>& theta,
                                         const std::vector<double>& x_r,
                                         const std::vector<int>& x_i,
                                         std::ostream* msgs = nullptr) {
    stan::math::check_not_nan("forward_pde", "theta", theta);
    const int need_sens = 0;
    std::vector<std::vector<double> > raw = pde(theta, need_sens, x_r, x_i, msgs);
    std::vector<double> res(raw.size());
    std::transform(raw.begin(), raw.end(),
                   res.begin(), [&theta](std::vector<double>& qoi_grad)
                   -> double { return qoi_grad[0]; } );
    return res;
  }

}
}

#endif
