#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_OBSERVER_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_OBSERVER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/arr/fun/sum.hpp>

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
template <typename F, typename T1, typename T2, typename T_t0, typename T_ts>
struct coupled_ode_observer {
  using return_t = return_type_t<T1, T2, T_t0, T_ts>;

  using ops_partials_t
      = operands_and_partials<std::vector<T1>, std::vector<T2>, T_t0, T_ts>;

  const F& f_;
  const std::vector<T1>& y0_;
  const T_t0& t0_;
  const std::vector<T_ts>& ts_;
  const std::vector<T2>& theta_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  std::ostream* msgs_;
  std::vector<std::vector<return_t>>& y_;
  const std::size_t N_;
  const std::size_t M_;
  const std::size_t index_offset_theta_;
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
  coupled_ode_observer(const F& f, const std::vector<T1>& y0,
                       const std::vector<T2>& theta, const T_t0& t0,
                       const std::vector<T_ts>& ts,
                       const std::vector<double>& x,
                       const std::vector<int>& x_int, std::ostream* msgs,
                       std::vector<std::vector<return_t>>& y)
      : f_(f),
        y0_(y0),
        t0_(t0),
        ts_(ts),
        theta_(theta),
        x_(x),
        x_int_(x_int),
        msgs_(msgs),
        y_(y),
        N_(y0.size()),
        M_(theta.size()),
        index_offset_theta_(is_constant_all<T1>::value ? 0 : N_ * N_),
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
  void operator()(const std::vector<double>& coupled_state, double t) {
    check_less("coupled_ode_observer", "time-state number", next_ts_index_,
               ts_.size());

    std::vector<return_t> yt;
    yt.reserve(N_);

    ops_partials_t ops_partials(y0_, theta_, t0_, ts_[next_ts_index_]);

    std::vector<double> dy_dt;
    if (!is_constant_all<T_ts>::value) {
      std::vector<double> y_dbl(coupled_state.begin(),
                                coupled_state.begin() + N_);
      dy_dt = f_(value_of(ts_[next_ts_index_]), y_dbl, value_of(theta_), x_,
                 x_int_, msgs_);
      check_size_match("coupled_ode_observer", "dy_dt", dy_dt.size(), "states",
                       N_);
    }

    for (size_t j = 0; j < N_; j++) {
      // iterate over parameters for each equation
      if (!is_constant_all<T1>::value) {
        for (std::size_t k = 0; k < N_; k++) {
          ops_partials.edge1_.partials_[k] = coupled_state[N_ + N_ * k + j];
        }
      }

      if (!is_constant_all<T2>::value) {
        for (std::size_t k = 0; k < M_; k++) {
          ops_partials.edge2_.partials_[k]
              = coupled_state[N_ + index_offset_theta_ + N_ * k + j];
        }
      }

      if (!is_constant_all<T_ts>::value) {
        ops_partials.edge4_.partials_[0] = dy_dt[j];
      }

      yt.emplace_back(ops_partials.build(coupled_state[j]));
    }

    y_.emplace_back(yt);
    next_ts_index_++;
  }
};

}  // namespace math

}  // namespace stan

#endif
