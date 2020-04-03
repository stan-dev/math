#ifndef STAN_MATH_REV_FUNCTOR_COUPLED_ODE_SYSTEM_HPP
#define STAN_MATH_REV_FUNCTOR_COUPLED_ODE_SYSTEM_HPP

#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/functor/coupled_ode_system.hpp>
#include <stan/math/rev/functor/cvodes_utils.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stdexcept>
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
template <typename T_initial, typename T_t0, typename T_ts, typename F, typename... Args>
struct coupled_ode_system_impl<false, T_initial, T_t0, T_ts, F, Args...> {
  using T_Return = return_type_t<T_initial, T_t0, T_ts, Args...>;

  const F& f_;
  const std::vector<T_initial>& y0_;
  std::tuple<const Args&...> args_tuple_;
  const size_t y0_vars_;
  const size_t args_vars_;
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
  coupled_ode_system_impl(const F& f, const std::vector<T_initial>& y0,
			  std::ostream* msgs,
			  const Args&... args)
      : f_(f),
        y0_(y0),
        args_tuple_(args...),
        y0_vars_(internal::count_vars(y0_)),
        args_vars_(internal::count_vars(args...)),
        N_(y0.size()),
        msgs_(msgs) {}

  /**
   * Calculates the derivative of the coupled ode system with respect
   * to time.
   *
   * This method uses nested autodiff and is not thread safe.
   *
   * @param[in] z state of the coupled ode system; this must be size
   *   <code>size()</code>
   * @param[out] dz_dt a vector of size <code>size()</code> with the
   *    derivatives of the coupled system with respect to time
   * @param[in] t time
   * @throw exception if the base ode function does not return the
   *    expected number of derivatives, N.
   */
  void operator()(const std::vector<double>& z, std::vector<double>& dz_dt,
                  double t) const {
    using std::vector;

    // Run nested autodiff in this scope
    nested_rev_autodiff nested;

    const vector<var> y_vars(z.begin(), z.begin() + N_);
    auto local_args_tuple = apply(
				  [&](auto&&... args) {
				    return std::tuple<decltype(internal::deep_copy(args))...>(
											      internal::deep_copy(args)...);
				  },
				  args_tuple_);
    
    vector<var> dy_dt_vars = apply(
				   [&](auto&&... args) { return f_(t, y_vars, msgs_, args...); },
				   local_args_tuple);

    check_size_match("coupled_ode_system", "dz_dt", dy_dt_vars.size(), "states",
                     N_);

    for (size_t i = 0; i < N_; i++) {
      dz_dt[i] = dy_dt_vars[i].val();
      dy_dt_vars[i].grad();
      for (size_t j = 0; j < y0_vars_; j++) {
	// orders derivatives by equation (i.e. if there are 2 eqns
	// (y1, y2) and 2 parameters (a, b), dy_dt will be ordered as:
	// dy1_dt, dy2_dt, dy1_da, dy2_da, dy1_db, dy2_db
	double temp_deriv = 0;
	const size_t offset = N_ + N_ * j;
	for (size_t k = 0; k < N_; k++) {
	  temp_deriv += z[N_ + N_ * j + k] * y_vars[k].adj();
	}
	
	dz_dt[N_ + N_ * j + i] = temp_deriv;
      }
      
      Eigen::VectorXd args_adjoints = Eigen::VectorXd::Zero(args_vars_);
      apply(
	    [&](auto&&... args) {
	      internal::accumulate_adjoints(args_adjoints.data(), args...);
	    },
	    local_args_tuple);
      for (size_t j = 0; j < args_vars_; j++) {
	double temp_deriv = args_adjoints(j);
	for (size_t k = 0; k < N_; k++) {
	  temp_deriv += z[N_ + N_ * y0_vars_ + N_ * j + k] * y_vars[k].adj();
	}
	
	dz_dt[N_ + N_ * y0_vars_ + N_ * j + i] = temp_deriv;
      }
      
      nested.set_zero_all_adjoints();
    }
  }

  template<typename T_t>
  std::vector<var> build_output(const std::vector<double>& dy0_dt0, const std::vector<double>& coupled_state, const T_t0& t0, const T_t& t) const {
    std::vector<double> y_dbl(coupled_state.data(), coupled_state.data() + N_);

    std::vector<var> yt;
    yt.reserve(N_);
    
    std::vector<double> dy_dt;
    if (is_var<T_t>::value) {
      std::vector<double> y_dbl(coupled_state.begin(),
                                coupled_state.begin() + N_);
      /*vector<var> dy_dt_vars = apply([&](const Args&... args) {
	  return f_(t, y_vars, args..., msgs_);
	  }, local_args_tuple);*/
      dy_dt = apply([&](auto&&... args) {
	  return f_(value_of(t), y_dbl, msgs_, value_of(args)...);
	}, args_tuple_);
      check_size_match("coupled_ode_observer", "dy_dt", dy_dt.size(), "states", N_);
    }
    
    for (size_t j = 0; j < N_; j++) {
      // When true this is 1 and not ts_.size() because there's
      //   only one time point involved with this output
      const size_t total_vars = y0_vars_ + args_vars_ + ((is_var<T_t>::value) ? 1 : 0) + ((is_var<T_t0>::value) ? 1 : 0);

      vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(total_vars);
      double* partials = ChainableStack::instance_->memalloc_.alloc_array<double>(total_vars);

      vari** varis_ptr = varis;
      double* partials_ptr = partials;

      // iterate over parameters for each equation
      varis_ptr = internal::save_varis(varis_ptr, y0_);
      for (std::size_t k = 0; k < y0_vars_; k++) {
	//*varis_ptr = y0[k].vi_;
	*partials_ptr = coupled_state[N_ + y0_vars_ * k + j];
	partials_ptr++;
      }

      varis_ptr = apply([&varis_ptr](auto&&... args) {
	  return internal::save_varis(varis_ptr, args...);
	}, args_tuple_);
      for (std::size_t k = 0; k < args_vars_; k++) {
	// dy[j]_dtheta[k]
	// theta[k].vi_
	*partials_ptr = coupled_state[N_ + N_ * y0_vars_ + N_ * k + j];
	partials_ptr++;
      }

      if ((stan::is_var<T_t>::value) ? 1 : 0) {
	varis_ptr = internal::save_varis(varis_ptr, t);
	*partials_ptr = dy_dt[j];
	partials_ptr++;
	// dy[j]_dcurrent_t
      }

      if ((stan::is_var<T_t0>::value) ? 1 : 0) {
	varis_ptr = internal::save_varis(varis_ptr, t0);
	*partials_ptr = -dy0_dt0[j];
	partials_ptr++;
	// dy[j]_dcurrent_t
      }

      yt.emplace_back(new precomputed_gradients_vari(coupled_state[j], total_vars, varis, partials));
    }

    return yt;
  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  size_t size() const { return N_ + N_ * y0_vars_ + N_ * args_vars_; }

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
  std::vector<double> initial_state() const {
    std::vector<double> initial(size(), 0.0);
    for (size_t i = 0; i < N_; i++) {
      initial[i] = value_of(y0_[i]);
    }
    for (size_t i = 0; i < y0_vars_; i++) {
      initial[N_ + i * N_ + i] = 1.0;
    }
    return initial;
  }
};

}  // namespace math
}  // namespace stan
#endif
