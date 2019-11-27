#ifndef STAN_MATH_REV_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_HPP
#define STAN_MATH_REV_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/core.hpp>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

/**
 * The <code>coupled_ode_system</code> template specialization
 * for known initial values and unknown parameters.
 *
 * <p>This coupled ode system has N + N * M states where N is the size of
 * the base ode system and M is the number of parameters.
 *
 * <p>For the coupled ode system, the first N states are the base
 * system's states: \f$ \frac{d x_n}{dt} \f$.
 *
 * <p>The next M states correspond to the sensitivities of the
 * parameters with respect to the first base system equation:
 * \f[
 *  \frac{d x_{N + m}}{dt}
 *  = \frac{d}{dt} \frac{\partial x_1}{\partial \theta_m}
 * \f]
 * for \f$ m \in {1, \ldots, M} \f$].
 *
 * <p> The next M states correspond to the sensitivites with respect
 * to the second base system equation, and so on through the last base
 * system equation.
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
 *   <code>operator()(double t, std::vector<double> y, std::vector<var> theta,
 *          std::vector<double> x, std::vector<int>x_int, std::ostream*
 * msgs)</code>
 */
template <typename F>
struct coupled_ode_system<F, double, var> {
  const F& f_;
  const std::vector<double>& y0_dbl_;
  const std::vector<var>& theta_;
  std::vector<var> theta_nochain_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  const size_t N_;
  const size_t M_;
  const size_t size_;
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
  coupled_ode_system(const F& f, const std::vector<double>& y0,
                     const std::vector<var>& theta,
                     const std::vector<double>& x,
                     const std::vector<int>& x_int, std::ostream* msgs)
      : f_(f),
        y0_dbl_(y0),
        theta_(theta),
        x_(x),
        x_int_(x_int),
        N_(y0.size()),
        M_(theta.size()),
        size_(N_ + N_ * M_),
        msgs_(msgs) {
    for (const var& p : theta) {
      theta_nochain_.emplace_back(var(new vari(p.val(), false)));
    }
  }

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

    try {
      start_nested();

      const vector<var> y_vars(z.begin(), z.begin() + N_);

      vector<var> dy_dt_vars = f_(t, y_vars, theta_nochain_, x_, x_int_, msgs_);

      check_size_match("coupled_ode_system", "dz_dt", dy_dt_vars.size(),
                       "states", N_);

      for (size_t i = 0; i < N_; i++) {
        dz_dt[i] = dy_dt_vars[i].val();
        dy_dt_vars[i].grad();

        for (size_t j = 0; j < M_; j++) {
          // orders derivatives by equation (i.e. if there are 2 eqns
          // (y1, y2) and 2 parameters (a, b), dy_dt will be ordered as:
          // dy1_dt, dy2_dt, dy1_da, dy2_da, dy1_db, dy2_db
          double temp_deriv = theta_nochain_[j].adj();
          const size_t offset = N_ + N_ * j;
          for (size_t k = 0; k < N_; k++) {
            temp_deriv += z[offset + k] * y_vars[k].adj();
          }

          dz_dt[offset + i] = temp_deriv;
        }

        set_zero_all_adjoints_nested();
        // Parameters stored on the outer (non-nested) nochain stack are not
        // reset to zero by the last call. This is done as a separate step here.
        // See efficiency note above on template specalization for more details
        // on this.
        for (size_t j = 0; j < M_; ++j) {
          theta_nochain_[j].vi_->set_zero_adjoint();
        }
      }
    } catch (const std::exception& e) {
      recover_memory_nested();
      throw;
    }
    recover_memory_nested();
  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  size_t size() const { return size_; }

  /**
   * Returns the initial state of the coupled system.
   *
   * The initial state of the coupled ode system is the same
   * as the base ode system. This is because the initial values
   * are known.
   *
   * There are N + N * M coupled states, where N is the number of base
   * ode system states and M is the number of parameters. The first N
   * correspond to the initial values provided. The next N * M states
   * are all 0.
   *
   * @return the initial condition of the coupled system, a vector of
   *   size N + N * M.
   */
  std::vector<double> initial_state() const {
    std::vector<double> state(size_, 0.0);
    for (size_t n = 0; n < N_; n++) {
      state[n] = y0_dbl_[n];
    }
    return state;
  }
};

/**
 * The <code>coupled_ode_system</code> template specialization
 * for unknown initial values and known parameters.
 *
 * <p>This coupled ode system has N + N * N states where N is the size
 * of the base ode system.
 *
 * <p>For the coupled ode system, the first N states are the base
 * system's states: \f$ \frac{d x_n}{dt} \f$.
 *
 * <p>The next N states correspond to the sensitivities of the
 * initial conditions with respect to the first base system equation:
 * \f[
 *  \frac{d x_{N + n}}{dt}
 *  = \frac{d}{dt} \frac{\partial x_1}{\partial y0_m}
 * \f]
 * for \f$ n \in {1, \ldots, N} \f$].
 *
 * <p> The next N states correspond to the sensitivites with respect
 * to the second base system equation, and so on through the last base
 * system equation.
 *
 * @tparam F base ode system functor. Must provide
 *   <code>operator()(double t, std::vector<var> y, std::vector<double> theta,
 *          std::vector<double> x, std::vector<int>x_int, std::ostream*
 * msgs)</code>
 */
template <typename F>
struct coupled_ode_system<F, var, double> {
  const F& f_;
  const std::vector<var>& y0_;
  const std::vector<double>& theta_dbl_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  std::ostream* msgs_;
  const size_t N_;
  const size_t M_;
  const size_t size_;

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
  coupled_ode_system(const F& f, const std::vector<var>& y0,
                     const std::vector<double>& theta,
                     const std::vector<double>& x,
                     const std::vector<int>& x_int, std::ostream* msgs)
      : f_(f),
        y0_(y0),
        theta_dbl_(theta),
        x_(x),
        x_int_(x_int),
        msgs_(msgs),
        N_(y0.size()),
        M_(theta.size()),
        size_(N_ + N_ * N_) {}

  /**
   * Calculates the derivative of the coupled ode system with respect
   * to time.
   *
   * This method uses nested autodiff and is not thread safe.
   *
   * @param[in] z state of the coupled ode syste; this must be
   *   size <code>size()</code>
   * @param[out] dz_dt a vector of length size() with the
   *   derivatives of the coupled system with respect to time
   * @param[in] t time
   * @throw exception if the base ode function does not return the
   *    expected number of derivatives, N.
   */
  void operator()(const std::vector<double>& z, std::vector<double>& dz_dt,
                  double t) const {
    using std::vector;

    try {
      start_nested();

      const vector<var> y_vars(z.begin(), z.begin() + N_);

      vector<var> dy_dt_vars = f_(t, y_vars, theta_dbl_, x_, x_int_, msgs_);

      check_size_match("coupled_ode_system", "dz_dt", dy_dt_vars.size(),
                       "states", N_);

      for (size_t i = 0; i < N_; i++) {
        dz_dt[i] = dy_dt_vars[i].val();
        dy_dt_vars[i].grad();

        for (size_t j = 0; j < N_; j++) {
          // orders derivatives by equation (i.e. if there are 2 eqns
          // (y1, y2) and 2 initial conditions (y0_a, y0_b), dy_dt will be
          // ordered as: dy1_dt, dy2_dt, dy1_d{y0_a}, dy2_d{y0_a}, dy1_d{y0_b},
          // dy2_d{y0_b}
          double temp_deriv = 0;
          const size_t offset = N_ + N_ * j;
          for (size_t k = 0; k < N_; k++) {
            temp_deriv += z[offset + k] * y_vars[k].adj();
          }

          dz_dt[offset + i] = temp_deriv;
        }

        set_zero_all_adjoints_nested();
      }
    } catch (const std::exception& e) {
      recover_memory_nested();
      throw;
    }
    recover_memory_nested();
  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  size_t size() const { return size_; }

  /**
   * Returns the initial state of the coupled system.
   *
   * <p>Because the starting state is unknown, the coupled system
   * incorporates the initial conditions as parameters.  At the initial
   * time the Jacobian of the integrated function is the identity matrix.
   *
   * @return the initial condition of the coupled system.
   *   This is a vector of length size() where the first N values are
   *   the initial condition of the base ODE and the remainder
   *   correspond to the identity matrix which is the Jacobian of the
   *   integrated function at the initial time-point.
   */
  std::vector<double> initial_state() const {
    std::vector<double> initial(size_, 0.0);
    for (size_t i = 0; i < N_; i++) {
      initial[i] = value_of(y0_[i]);
    }
    for (size_t i = 0; i < N_; i++) {
      initial[N_ + i * N_ + i] = 1.0;
    }
    return initial;
  }
};

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
template <typename F>
struct coupled_ode_system<F, var, var> {
  const F& f_;
  const std::vector<var>& y0_;
  const std::vector<var>& theta_;
  std::vector<var> theta_nochain_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  const size_t N_;
  const size_t M_;
  const size_t size_;
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
  coupled_ode_system(const F& f, const std::vector<var>& y0,
                     const std::vector<var>& theta,
                     const std::vector<double>& x,
                     const std::vector<int>& x_int, std::ostream* msgs)
      : f_(f),
        y0_(y0),
        theta_(theta),
        x_(x),
        x_int_(x_int),
        N_(y0.size()),
        M_(theta.size()),
        size_(N_ + N_ * (N_ + M_)),
        msgs_(msgs) {
    for (const var& p : theta) {
      theta_nochain_.emplace_back(var(new vari(p.val(), false)));
    }
  }

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

    try {
      start_nested();

      const vector<var> y_vars(z.begin(), z.begin() + N_);

      vector<var> dy_dt_vars = f_(t, y_vars, theta_nochain_, x_, x_int_, msgs_);

      check_size_match("coupled_ode_system", "dz_dt", dy_dt_vars.size(),
                       "states", N_);

      for (size_t i = 0; i < N_; i++) {
        dz_dt[i] = dy_dt_vars[i].val();
        dy_dt_vars[i].grad();

        for (size_t j = 0; j < N_; j++) {
          // orders derivatives by equation (i.e. if there are 2 eqns
          // (y1, y2) and 2 parameters (a, b), dy_dt will be ordered as:
          // dy1_dt, dy2_dt, dy1_da, dy2_da, dy1_db, dy2_db
          double temp_deriv = 0;
          const size_t offset = N_ + N_ * j;
          for (size_t k = 0; k < N_; k++) {
            temp_deriv += z[offset + k] * y_vars[k].adj();
          }

          dz_dt[offset + i] = temp_deriv;
        }

        for (size_t j = 0; j < M_; j++) {
          double temp_deriv = theta_nochain_[j].adj();
          const size_t offset = N_ + N_ * N_ + N_ * j;
          for (size_t k = 0; k < N_; k++) {
            temp_deriv += z[offset + k] * y_vars[k].adj();
          }

          dz_dt[offset + i] = temp_deriv;
        }

        set_zero_all_adjoints_nested();
        // Parameters stored on the outer (non-nested) nochain stack are not
        // reset to zero by the last call. This is done as a separate step here.
        // See efficiency note above on template specalization for more details
        // on this.
        for (size_t j = 0; j < M_; ++j) {
          theta_nochain_[j].vi_->set_zero_adjoint();
        }
      }
    } catch (const std::exception& e) {
      recover_memory_nested();
      throw;
    }
    recover_memory_nested();
  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  size_t size() const { return size_; }

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
    std::vector<double> initial(size_, 0.0);
    for (size_t i = 0; i < N_; i++) {
      initial[i] = value_of(y0_[i]);
    }
    for (size_t i = 0; i < N_; i++) {
      initial[N_ + i * N_ + i] = 1.0;
    }
    return initial;
  }
};

}  // namespace math
}  // namespace stan
#endif
