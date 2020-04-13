#ifndef STAN_MATH_PRIM_FUNCTOR_COUPLED_ODE_SYSTEM_HPP
#define STAN_MATH_PRIM_FUNCTOR_COUPLED_ODE_SYSTEM_HPP

#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * The <code>coupled_ode_system</code> template specialization
 * for unknown initial values and unknown parameters.
 *
 * <p>This coupled ode system has N + (N +  M) * N states where N is
 * the size of the base ode system and M is the number of parameters.
 *
 * <p>For the coupled ode system, the first N states are the base
 * system's states: \f$ \frac{d x_n}{dt} \f$.
 *
 * <p>The next N + M states correspond to the sensitivities of the
 * initial conditions, then to the sensitivities of the parameters
 * with respect to the to the first base system equation:
 *
 * \f[
 *   \frac{d x_{N + n}}{dt}
 *     = \frac{d}{dt} \frac{\partial x_1}{\partial y0_n}
 * \f]
 *
 * \f[
 *   \frac{d x_{N + N + m}}{dt}
 *     = \frac{d}{dt} \frac{\partial x_1}{\partial \theta_m}
 * \f]
 *
 * <p>The next N + M states correspond to the sensitivities
 * of the initial conditions followed by the sensitivites of the
 * parameters with respect to the second base system equation, and
 * so on through the last base system equation.
 *
 * <p>Note: Calculating the sensitivity system requires the Jacobian
 * of the base ODE RHS wrt to the parameters theta. The parameter
 * vector theta is constant for successive calls to the exposed
 * operator(). For this reason, the parameter vector theta is copied
 * upon construction onto the nochain var autodiff tape which is used
 * in the the nested autodiff performed in the operator() of this
 * adaptor. Doing so reduces the size of the nested autodiff and
 * speeds up autodiff. As a side effect, the parameter vector theta
 * will remain on the nochain autodiff part of the autodiff tape being
 * in use even after destruction of the given instance. Moreover, the
 * adjoint zeroing for the nested system does not cover the theta
 * parameter vector part of the nochain autodiff tape and is therefore
 * set to zero using a dedicated loop.
 *
 * @tparam F base ode system functor. Must provide
 *   <code>operator()(double t, std::vector<var> y, std::vector<var> theta,
 *          std::vector<double> x, std::vector<int>x_int, std::ostream*
 * msgs)</code>
 */

template <bool Enable, typename F, typename T_initial, typename... Args>
struct coupled_ode_system_impl;
  
template <typename F, typename T_initial, typename... Args>
struct coupled_ode_system_impl<true, F, T_initial, Args...> {
  const F& f_;
  const std::vector<double>& y0_;
  std::tuple<const Args&...> args_tuple_;
  const size_t N_;
  std::ostream* msgs_;

  /**
   * Construct a coupled ode system from the base system function,
   * initial state of the base system, parameters, and a stream for
   * messages.
   *
   * @param[in] f the base ODE system functor
   * @param[in] y0 the initial state of the base ode
   * @param[in] theta parameters of the base ode
   * @param[in] x real data
   * @param[in] x_int integer data
   * @param[in, out] msgs stream for messages
   */
  coupled_ode_system_impl(const F& f, const std::vector<double>& y0,
			  std::ostream* msgs, const Args&... args)
      : f_(f),
        y0_(y0),
	args_tuple_(args...),
        N_(y0.size()),
        msgs_(msgs) {
  }

  void operator()(const std::vector<double>& y, std::vector<double>& dy_dt,
                  double t) const {
    dy_dt = apply([&](const Args&... args) {
	return f_(t, y, msgs_, args...);
      }, args_tuple_);

    check_size_match("coupled_ode_system", "y", y.size(), "dy_dt",
                     dy_dt.size());
  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  size_t size() const { return N_; }

  /**
   * Returns the initial state of the coupled system.
   *
   * <p>Because the starting state is unknown, the coupled system
   * incorporates the initial conditions as parameters. At the initial
   * time the Jacobian of the integrated function is the identity
   * matrix. In addition the coupled system includes the Jacobian of
   * the integrated function wrt to the parameters. This Jacobian is
   * zero at the initial time-point.
   *
   * @return the initial condition of the coupled system.  This is a
   *   vector of length size() where the first N values are the
   *   initial condition of the base ODE and the next N*N elements
   *   correspond to the identity matrix which is the Jacobian of the
   *   integrated function at the initial time-point. The last N*M
   *   elements are all zero as these are the Jacobian wrt to the
   *   parameters at the initial time-point, which is zero.
   */
  std::vector<double> initial_state() const { return value_of(y0_); }
};

template <typename F, typename T_initial, typename... Args>
struct coupled_ode_system :
    public coupled_ode_system_impl<std::is_arithmetic<return_type_t<T_initial, Args...>>::value,
				   F, T_initial, Args...> {

  coupled_ode_system(const F& f, const std::vector<T_initial>& y0,
		     std::ostream* msgs, const Args&... args)
    : coupled_ode_system_impl<std::is_arithmetic<return_type_t<T_initial, Args...>>::value,
			      F, T_initial, Args...>(f, y0, msgs, args...) {}
};

}  // namespace math
}  // namespace stan
#endif
