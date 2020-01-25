#ifndef STAN_MATH_REV_FUNCTOR_COUPLED_ODE_OBSERVER_HPP
#define STAN_MATH_REV_FUNCTOR_COUPLED_ODE_OBSERVER_HPP

#include <stan/math/rev/functor/cvodes_utils.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/sum.hpp>

#include <vector>

namespace stan {
namespace math {

/**
 * Observer for the coupled states.  Holds a reference to an
 * externally defined vector of vectors passed in at construction time
 * which holds the final result on the AD stack. Thus, whenever any of
 * the inputs is varying, then the output will be varying as well. The
 * sensitivities of the initials and the parameters are taken from the
 * coupled state in the order as defined by the
 * coupled_ode_system. The sensitivities at each time-point is simply
 * the ODE RHS evaluated at that time point.
 *
 * The output of this class is for all time-points in the ts vector
 * which does not contain the initial time-point by the convention
 * used in stan-math.
 *
 */
template <typename T1, typename T_t0, typename T_ts, typename F, typename... Args>
struct coupled_ode_observer {
  using ReturnType = return_type_t<T1, T_t0, T_ts, Args...>;
  
  const F& f_;
  const std::vector<T1>& y0_;
  const T_t0& t0_;
  const std::vector<T_ts>& ts_;
  std::ostream* msgs_;
  std::vector<std::vector<ReturnType>>& y_;
  const std::size_t N_;
  std::tuple<const Args&...> args_tuple_;
  int next_ts_index_;

  /**
   * Construct a coupled ODE observer for the specified coupled
   * vector.
   *
   * @tparam F type of ODE system function.
   * @tparam T1 type of scalars for initial values.
   * @tparam T2 type of scalars for parameters.
   * @tparam T_t0 type of scalar of initial time point.
   * @tparam T_ts type of time-points where ODE solution is returned.
   * @param[in] f functor for the base ordinary differential equation.
   * @param[in] y0 initial state.
   * @param[in] theta parameter vector for the ODE.
   * @param[in] t0 initial time.
   * @param[in] ts times of the desired solutions, in strictly
   * increasing order, all greater than the initial time.
   * @param[in] x continuous data vector for the ODE.
   * @param[in] x_int integer data vector for the ODE.
   * @param[out] msgs the print stream for warning messages.
   * @param[out] y reference to a vector of vector of the final return
   */
  coupled_ode_observer(const F& f,
		       const std::vector<T1>& y0,
		       const T_t0& t0,
                       const std::vector<T_ts>& ts,
		       const Args&... args,
		       std::ostream* msgs,
                       std::vector<std::vector<ReturnType>>& y)
      : f_(f),
        y0_(y0),
        t0_(t0),
        ts_(ts),
        msgs_(msgs),
        y_(y),
	args_tuple_(args...),
        N_(y0.size()),
        next_ts_index_(0) {}

  /**
   * Callback function for ODE solvers to record values. The coupled
   * state returned from the solver is added directly to the AD tree.
   *
   * The coupled state follows the convention as defined in the
   * coupled_ode_system. In brief, the coupled state consists of {f,
   * df/dy0, df/dtheta}. Here df/dy0 and df/dtheta are only present if
   * their respective sensitivites have been requested.
   *
   * @param coupled_state solution at the specified time.
   * @param t time of solution. The time must correspond to the
   * element ts[next_ts_index_]
   */
  template<bool return_type_is_var = is_var<ReturnType>::value>
  std::enable_if_t<return_type_is_var>
  operator()(const std::vector<double>& coupled_state, double t) {
    check_less("coupled_ode_observer", "time-state number", next_ts_index_,
               ts_.size());
    
    std::vector<ReturnType> yt;
    yt.reserve(N_);

    std::vector<double> dy_dt;
    if (stan::is_var<T_ts>::value) {
      std::vector<double> y_dbl(coupled_state.begin(),
                                coupled_state.begin() + N_);
      dy_dt = apply([&](const Args&... args) {
	  return f_(value_of(ts_[next_ts_index_]), y_dbl, value_of(args)..., msgs_);
	}, args_tuple_);
      check_size_match("coupled_ode_observer", "dy_dt", dy_dt.size(), "states", N_);
    }

    const size_t y0_vars = internal::count_vars(y0_);
    const size_t args_vars = apply([&](const Args&... args) { return internal::count_vars(args...); }, args_tuple_);
    // When true this is 1 and not ts_.size() because there's
    //   only one time point involved with this output
    const size_t time_vars = (stan::is_var<T_ts>::value) ? 1 : 0;
    const size_t total_vars = y0_vars + args_vars + time_vars;
    
    for (size_t j = 0; j < N_; j++) {
      vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(total_vars);
      double* partials = ChainableStack::instance_->memalloc_.alloc_array<double>(total_vars);

      vari** varis_ptr = varis;
      double* partials_ptr = partials;

      size_t offset = 0;

      // iterate over parameters for each equation
      varis_ptr = internal::save_varis(varis_ptr, y0_);
      for (std::size_t k = 0; k < y0_vars; k++) {
	//*varis_ptr = y0[k].vi_;
	*partials_ptr = coupled_state[N_ + y0_vars * k + j];
	partials_ptr++;
      }

      varis_ptr = apply([&varis_ptr](const Args&... args) {
	  return internal::save_varis(varis_ptr, args...);
	}, args_tuple_);
      for (std::size_t k = 0; k < args_vars; k++) {
	// dy[j]_dtheta[k]
	// theta[k].vi_
	*partials_ptr = coupled_state[N_ + N_ * y0_vars + N_ * k + j];
	partials_ptr++;
      }

      if (time_vars > 0) {
	varis_ptr = internal::save_varis(varis_ptr, ts_[next_ts_index_]);
	*partials_ptr = dy_dt[j];
	partials_ptr++;
	// dy[j]_dcurrent_t
      }

      yt.emplace_back(new precomputed_gradients_vari(coupled_state[j], total_vars, varis, partials));
    }

    y_.emplace_back(yt);
    next_ts_index_++;
  }

  template<bool return_type_is_var = is_var<ReturnType>::value>
  std::enable_if_t<not return_type_is_var>
  operator()(const std::vector<double>& coupled_state, double t) {
    y_.emplace_back(coupled_state);
  }
};

}  // namespace math

}  // namespace stan

#endif
