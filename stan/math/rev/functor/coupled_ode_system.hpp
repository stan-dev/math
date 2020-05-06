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
template <typename F, typename T_initial, typename... Args>
struct coupled_ode_system_impl<false, F, T_initial, Args...> {
  const F& f_;
  const Eigen::Matrix<T_initial, Eigen::Dynamic, 1>& y0_;
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
  coupled_ode_system_impl(const F& f, const Eigen::Matrix<T_initial, Eigen::Dynamic, 1>& y0,
                          std::ostream* msgs, const Args&... args)
      : f_(f),
        y0_(y0),
        args_tuple_(args...),
        y0_vars_(count_vars(y0_)),
        args_vars_(count_vars(args...)),
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
  void operator()(const Eigen::VectorXd& z, Eigen::VectorXd& dz_dt,
                  double t) const {
    using std::vector;

    dz_dt.resize(size());
    
    // Run nested autodiff in this scope
    nested_rev_autodiff nested;

    Eigen::Matrix<var, Eigen::Dynamic, 1> y_vars(N_);
    for(size_t n = 0; n < N_; ++n)
      y_vars(n) = z(n);

    auto local_args_tuple = apply(
        [&](auto&&... args) {
          return std::tuple<decltype(deep_copy_vars(args))...>(
              deep_copy_vars(args)...);
        },
        args_tuple_);

    Eigen::Matrix<var, Eigen::Dynamic, 1> f_y_t_vars
      = apply([&](auto&&... args) { return f_(t, y_vars, msgs_, args...); },
	      local_args_tuple);

    check_size_match("coupled_ode_system", "dy_dt", f_y_t_vars.size(), "states",
                     N_);

    for (size_t i = 0; i < N_; i++) {
      dz_dt(i) = f_y_t_vars(i).val();
      f_y_t_vars(i).grad();
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
            accumulate_adjoints(args_adjoints.data(), args...);
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
  Eigen::VectorXd initial_state() const {
    Eigen::VectorXd initial = Eigen::VectorXd::Zero(size());
    for (size_t i = 0; i < N_; i++) {
      initial(i) = value_of(y0_[i]);
    }
    for (size_t i = 0; i < y0_vars_; i++) {
      initial(N_ + i * N_ + i) = 1.0;
    }
    return initial;
  }
};

}  // namespace math
}  // namespace stan
#endif
