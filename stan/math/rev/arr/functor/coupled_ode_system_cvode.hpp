#ifndef STAN_MATH_REV_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_CVODE_HPP
#define STAN_MATH_REV_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_CVODE_HPP

#include <stan/math/rev/arr/functor/coupled_ode_system.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system_cvode.hpp>

#include <cvode/cvode_band.h>

namespace stan {
  namespace math {

    /**
     * The coupled ODE system for known initial values and unknown
     * parameters.
     *
     * <p>If the base ODE state is size N and there are M parameters,
     * the coupled system has N + N * M states.
     * <p>The first N states correspond to the base system's N states:
     * \f$ \frac{d x_n}{dt} \f$
     *
     * <p>The next M states correspond to the sensitivities of the
     * parameters with respect to the first base system equation:
     * \f[
     *   \frac{d x_{N+m}}{dt}
     *   = \frac{d}{dt} \frac{\partial x_1}{\partial \theta_m}
     * \f]
     *
     * <p>The final M states correspond to the sensitivities with respect
     * to the second base system equation, etc.
     *
     * @tparam F type of functor for the base ode system.
     */
    template <typename F>
    class coupled_ode_system_cvode <F, double, stan::math::var>:
      public coupled_ode_system<F, double, stan::math::var> {
      /**
       * Construct a coupled ODE system with the specified base
       * ODE system, base initial state, parameters, data, and a
       * message stream.
       *
       * @param[in] f the base ODE system functor.
       * @param[in] y0 the initial state of the base ode.
       * @param[in] theta parameters of the base ode.
       * @param[in] x real data.
       * @param[in] x_int integer data.
       * @param[in, out] msgs stream to which messages are printed.
       */
      coupled_ode_system_cvode(const F& f,
                               const std::vector<double>& y0,
                               const std::vector<stan::math::var>& theta,
                               const std::vector<double>& x,
                               const std::vector<int>& x_int,
                               double rel_tol,
                               double abs_tol,
                               long int max_num_steps,
                               std::ostream* msgs)
        : coupled_ode_system(f, y0, theat, x, x_int, msg),
          cvode_mem_(NULL),
          state_(N_ * (M_ + 1),
          cvode_state_(N_ * (M_ + 1), &state[0]) {
            instantiate_cvode_mem(rel_tol, abs_tol, max_num_steps);
            check_flag(CVBand(cvode_mem_, N_ * (M_ + 1), N_ - 1, N_ - 1),
                       "CVBand");
            check_flag(CVDlsSetBandJacFn(cvode_mem_, &ode::banded_jacobian),
                       "CVDlsSetBandJacFn");
      }

      // Forward the CVode array-based call from
      // the ODE RHS to the Stan vector-based call
      void operator()(const double y[], double dy_dt[], double t) {
        std::vector<double> y_vec(N_ * (M_ + 1));
        std::copy(y, y + N_, y_vec);

        std::vector<double> dy_dt_vec(N_ * (M_ + 1));
        (*this)(y_vec, dy_dt_vec, t);

        std::copy(dy_dt_vec.begin(), dy_dt_vec.end(), dy_dt);
      }

      // Static wrapper for CVode callback
      static int banded_jacobian(long int N, long int m_upper, long int m_lower,
                                realtype t, N_Vector y, N_Vector fy,
                                DlsMat J, void *J_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        ode* explicit_ode = static_cast<ode*>(J_data);
        return explicit_ode->banded_jacobian(NV_DATA_S(y), J->cols,
                                             long int s_mu, t);
      }

      // J[j][i] = d(ydot[i])/d(y[j])
      void banded_jacobian(const double y[], double *J[],
                          long int s_mu, double t) {
        using stan::math::var;

        std::vector<double> y_vec(N_);
        std::copy(y, y + N_, y_vec);

        std::vector<double> grad(N_);

        // Column major storage assumed for base Jacobian
        std::vector<std::vector<double> > base_J;

        try {
          stan::math::start_nested();

          std::vector<var> y_vars(y.begin(), y.end());
          std::vector<var> dy_dt_vars = f_(t, y_vars, theta_dbl_,
                                           x_, x_int_, msgs_);

          for (size_t i = 0; i < N_; i++) {
            set_zero_all_adjoints_nested();
            dy_dt_vars[i].grad(y_vars, grad);
            for (size_t j = 0; j < N_; ++j)
              base_J[j][i] = grad[j];
          }
        } catch (const std::exception& e) {
          stan::math::recover_memory_nested();
          throw;
        }
        stan::math::recover_memory_nested();
      }

      for (int j = 0; j < N_ * (M_ + 1); ++j) {
        std::memcpy(J[j] + s_mu - j % N_,
                    &(base_J[j % N_])(0),
                    N_ * sizeof(double));
      }
    };

    /**
     * The coupled ODE system for unknown initial values and known
     * parameters.
     *
     * <p>If the original ODE has states of size N, the
     * coupled system has N + N * N states. (derivatives of each
     * state with respect to each initial value)
     *
     * <p>The coupled system has N + N * N states, where N is the size of
     * the state vector in the base system.
     *
     * <p>The first N states correspond to the base system's N states:
     * \f$ \frac{d x_n}{dt} \f$
     *
     * <p>The next N states correspond to the sensitivities of the initial
     * conditions with respect to the to the first base system equation:
     * \f[
     *  \frac{d x_{N+n}}{dt}
     *     = \frac{d}{dt} \frac{\partial x_1}{\partial y0_n}
     * \f]
     *
     * <p>The next N states correspond to the sensitivities with respect
     * to the second base system equation, etc.
     *
     * @tparam F type of base ODE system functor
     */
    template <typename F>
    class coupled_ode_system_cvode <F, stan::math::var, double>:
      public coupled_ode_system<F, stan::math::var, double> {
      /**
       * Construct a coupled ODE system for an unknown initial state
       * and known parameters givne the specified base system functor,
       * base initial state, parameters, data, and an output stream
       * for messages.
       *
       * @param[in] f base ODE system functor.
       * @param[in] y0 initial state of the base ODE.
       * @param[in] theta system parameters.
       * @param[in] x real data.
       * @param[in] x_int integer data.
       * @param[in, out] msgs output stream for messages.
       */
      coupled_ode_system_cvode(const F& f,
                         const std::vector<stan::math::var>& y0,
                         const std::vector<double>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         double rel_tol,
                         double abs_tol,
                         long int max_num_steps,
                         std::ostream* msgs)
        : coupled_ode_system(f, y0, theat, x, x_int, msg),
          cvode_mem_(NULL),
          state_(N_ * (N_ + 1),
          cvode_state_(N_ * (N_ + 1), &state[0]) {
            instantiate_cvode_mem(rel_tol, abs_tol, max_num_steps);
            check_flag(CVBand(cvode_mem_, N_ * (N_ + 1), N_ - 1, N_ - 1),
                       "CVBand");
            check_flag(CVDlsSetBandJacFn(cvode_mem_, &ode::banded_jacobian),
                       "CVDlsSetBandJacFn");
      }

      // Forward the CVode array-based call from
      // the ODE RHS to the Stan vector-based call
      void operator()(const double y[], double dy_dt[], double t) {
        std::vector<double> y_vec(N_ * (N_ + 1));
        std::copy(y, y + N_, y_vec);

        std::vector<double> dy_dt_vec(N_ * (N_ + 1));
        (*this)(y_vec, dy_dt_vec, t);

        std::copy(dy_dt_vec.begin(), dy_dt_vec.end(), dy_dt);
      }

      // Static wrapper for CVode callback
      static int banded_jacobian(long int N, long int m_upper, long int m_lower,
                                realtype t, N_Vector y, N_Vector fy,
                                DlsMat J, void *J_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        ode* explicit_ode = static_cast<ode*>(J_data);
        return explicit_ode->banded_jacobian(NV_DATA_S(y), J->cols,
                                             long int s_mu, t);
      }

      // J[j][i] = d(ydot[i])/d(y[j])
      void banded_jacobian(const double y[], double *J[],
                          long int s_mu, double t) {
        using stan::math::var;

        std::vector<double> y_vec(N_);
        for (size_t n = 0; n < N_; n++)
          y_vec[n] = y[n] + y0_dbl_[n];

        std::vector<double> grad(N_);

        // Column major storage assumed for base Jacobian
        std::vector<std::vector<double> > base_J;

        try {
          stan::math::start_nested();

          std::vector<var> y_vars(y.begin(), y.end());
          std::vector<var> dy_dt_vars = f_(t, y_vars, theta_dbl_,
                                           x_, x_int_, msgs_);

          for (size_t i = 0; i < N_; i++) {
            set_zero_all_adjoints_nested();
            dy_dt_vars[i].grad(y_vars, grad);
            for (size_t j = 0; j < N_; ++j)
              base_J[j][i] = grad[j];
          }
        } catch (const std::exception& e) {
          stan::math::recover_memory_nested();
          throw;
        }
        stan::math::recover_memory_nested();
      }

      for (int j = 0; j < N_ * (N_ + 1); ++j) {
        std::memcpy(J[j] + s_mu - j % N_,
                    &(base_J[j % N_])(0),
                    N_ * sizeof(double));
      }
    };

    /**
     * The coupled ode system for unknown intial values and unknown
     * parameters.
     *
     * <p>The coupled system has N + N * (N + M) states, where N is
     * size of the base ODE state vector and M is the number of
     * parameters.
     *
     * <p>The first N states correspond to the base system's N states:
     *   \f$ \frac{d x_n}{dt} \f$
     *
     * <p>The next N+M states correspond to the sensitivities of the
     * initial conditions, then to the parameters with respect to the
     * to the first base system equation:
     *
     * \f[
     *   \frac{d x_{N + n}}{dt}
     *     = \frac{d}{dt} \frac{\partial x_1}{\partial y0_n}
     * \f]
     *
     * \f[
     *   \frac{d x_{N+N+m}}{dt}
     *     = \frac{d}{dt} \frac{\partial x_1}{\partial \theta_m}
     * \f]
     *
     * <p>The next N+M states correspond to the sensitivities with
     * respect to the second base system equation, etc.
     *
     * <p>If the original ode has a state vector of size N states and
     * a parameter vector of size M, the coupled system has N + N * (N
     * + M) states. (derivatives of each state with respect to each
     * initial value and each theta)
     *
     * @tparam F the functor for the base ode system
     */
    template <typename F>
    struct coupled_ode_system_cvode <F, stan::math::var, stan::math::var>:
      public coupled_ode_system <F, stan::math::var, stan::math::var> {
      /**
       * Construct a coupled ODE system with unknown initial value and
       * known parameters, given the base ODE system functor, the
       * initial state of the base ODE, the parameters, data, and an
       * output stream to which to write messages.
       *
       * @param[in] f the base ode system functor.
       * @param[in] y0 the initial state of the base ode.
       * @param[in] theta parameters of the base ode.
       * @param[in] x real data.
       * @param[in] x_int integer data.
       * @param[in, out] msgs output stream to which to print messages.
       */
      coupled_ode_system_cvode(const F& f,
                         const std::vector<stan::math::var>& y0,
                         const std::vector<stan::math::var>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         double rel_tol,
                         double abs_tol,
                         long int max_num_steps,
                         std::ostream* msgs)
        : coupled_ode_system(f, y0, theat, x, x_int, msg),
          cvode_mem_(NULL),
          state_(N_ * (M_ + N_ + 1),
          cvode_state_(N_ * (M_ + N_ + 1), &state[0]) {
            instantiate_cvode_mem(rel_tol, abs_tol, max_num_steps);
            check_flag(CVBand(cvode_mem_, N_ * (M_ + N_ + 1), N_ - 1, N_ - 1),
                       "CVBand");
            check_flag(CVDlsSetBandJacFn(cvode_mem_, &ode::banded_jacobian),
                       "CVDlsSetBandJacFn");
      }

      // Forward the CVode array-based call from
      // the ODE RHS to the Stan vector-based call
      void operator()(const double y[], double dy_dt[], double t) {
        std::vector<double> y_vec(N_ * (M_ + N_ + 1));
        std::copy(y, y + N_, y_vec);

        std::vector<double> dy_dt_vec(N_ * (M_ + N_ + 1));
        (*this)(y_vec, dy_dt_vec, t);

        std::copy(dy_dt_vec.begin(), dy_dt_vec.end(), dy_dt);
      }

      // Static wrapper for CVode callback
      static int banded_jacobian(long int N, long int m_upper, long int m_lower,
                                realtype t, N_Vector y, N_Vector fy,
                                DlsMat J, void *J_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        ode* explicit_ode = static_cast<ode*>(J_data);
        return explicit_ode->banded_jacobian(NV_DATA_S(y), J->cols,
                                             long int s_mu, t);
      }

      // J[j][i] = d(ydot[i])/d(y[j])
      void banded_jacobian(const double y[], double *J[],
                          long int s_mu, double t) {
        using stan::math::var;

        std::vector<double> y_vec(N_);
        for (size_t n = 0; n < N_; n++)
          y_vec[n] = y[n] + y0_dbl_[n];

        std::vector<double> grad(N_);

        // Column major storage assumed for base Jacobian
        std::vector<std::vector<double> > base_J;

        try {
          stan::math::start_nested();

          std::vector<var> y_vars(y.begin(), y.end());
          std::vector<var> dy_dt_vars = f_(t, y_vars, theta_dbl_,
                                           x_, x_int_, msgs_);

          for (size_t i = 0; i < N_; i++) {
            set_zero_all_adjoints_nested();
            dy_dt_vars[i].grad(y_vars, grad);
            for (size_t j = 0; j < N_; ++j)
              base_J[j][i] = grad[j];
          }
        } catch (const std::exception& e) {
          stan::math::recover_memory_nested();
          throw;
        }
        stan::math::recover_memory_nested();
      }

      for (int j = 0; j < N_ * (M_ + N_ + 1); ++j) {
        std::memcpy(J[j] + s_mu - j % N_,
                    &(base_J[j % N_])(0),
                    N_ * sizeof(double));
      }

    };
  } // math
} // stan

#endif
