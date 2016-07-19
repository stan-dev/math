#ifndef STAN_MATH_REV_MAT_FUNCTOR_COUPLED_ODE_SYSTEM_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_COUPLED_ODE_SYSTEM_HPP

#include <stan/math/rev/mat/functor/ode_system.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/functor/coupled_ode_system.hpp>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace stan {
  namespace math {

    /**
     * The coupled ode system suitable for varying/fixed
     * initials/parameters.
     *
     * <p> The base system has N states and M parameters.
     *
     * If any sensitivities are requested then the state vector gets
     * enlarged accordingly by the number of sensitivities requested
     * (for each sensitivity all N states are being calculated).
     *
     * If sensitivities for the initials are requested, then the state
     * vector is enlarged by N * N. The initial state for the
     * sensitivities of the initials is the identity matrix.
     *
     * In case the M parameter sensitivities are requested, the state
     * vector is enhanced by N * M states which are initialized at t0
     * to 0.
     *
     * <p>The largest possible coupled system has N + N * (N + M)
     * states, where N is size of the base ODE state vector and M is
     * the number of parameters.
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
     * + M) states (derivatives of each state with respect to each
     * initial value and each theta).
     *
     * @tparam F the functor for the base ode system
     */
    template <typename F, typename T_initial, typename T_param>
    class coupled_ode_system {
      const std::vector<T_initial>& y0_;
      const std::vector<T_param>& theta_;
      const size_t N_;
      const size_t M_;
      const size_t S_;
      const size_t size_;
      const ode_system<F> ode_system_;

      typedef stan::is_var<T_initial> initial_var;
      typedef stan::is_var<T_param> param_var;
      typedef typename stan::return_type<T_initial, T_param>::type return_t;

    public:
      /**
       * Construct a coupled ODE system with fixed/varying
       * initials/parameters, given the base ODE system functor, the
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
      coupled_ode_system(const F& f,
                         const std::vector<T_initial>& y0,
                         const std::vector<T_param>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         std::ostream* msgs)
        : y0_(y0),
          theta_(theta),
          N_(y0.size()),
          M_(theta.size()),
          S_((initial_var::value ? N_ : 0) + (param_var::value ? M_ : 0)),
          size_(N_ * (1 + S_)),
          ode_system_(f, stan::math::value_of(theta), x, x_int, msgs) { }

      /**
       * Populates the derivative vector with derivatives of the
       * coupled ODE system state with respect to time evaluated at the
       * specified state and specified time.
       *
       * @param[in]  z the current state of the coupled ode system,
       * of size <code>size()</code>.
       * @param[in, out] dz_dt populate with the derivatives of the
       * coupled system evaluated at the specified state and time.
       * @param[in] t time.
       * @throw exception if the base system does not return a
       * derivative vector of the same size as the state vector.
       *
       * y is the base ODE system state
       *
       */
      inline void
      operator()(const std::vector<double>& z,
                 std::vector<double>& dz_dt,
                 double t) const {
        rhs_sens(y0_, theta_, z, dz_dt, t);
      }

      /**
       * Returns the size of the coupled system.
       *
       * @return size of the coupled system.
       */
      size_t size() const {
        return size_;
      }

      /**
       * Returns the initial state of the coupled system.
       *
       * In case the initials are varying, these get assigned the
       * identity matrix as initial.
       *
       * @return the initial condition of the coupled system.  This is
       * a vector of length size().
       */
      std::vector<double> initial_state() const {
        std::vector<double> initial(size_, 0.0);
        for (size_t i = 0; i < N_; i++)
          initial[i] = stan::math::value_of(y0_[i]);
        if (initial_var::value) {
          for (size_t i = 0; i < N_; i++)
            initial[N_ + i * N_ + i] = 1.0;
        }
        return initial;
      }

    private:
      inline void
      rhs_sens(const std::vector<stan::math::var>& initial,
               const std::vector<stan::math::var>& param,
               const std::vector<double>& z,
               std::vector<double>& dz_dt,
               double t) const {
        using Eigen::Map;
        using Eigen::MatrixXd;
        using Eigen::VectorXd;

        const std::vector<double> y_base(z.begin(), z.begin() + N_);
        Map<VectorXd> dz_dt_eig(&dz_dt[0], N_);
        Map<MatrixXd> dZ_dt_sens(&dz_dt[0] + N_, N_, S_);
        Map<const MatrixXd> Z_sens(&z[0] + N_, N_, S_);
        // write Jtheta directly into correct position of dZ_dt
        Map<MatrixXd> dZ_dt_sens_param(&dz_dt[0] + N_ + N_ * N_, N_, M_);
        MatrixXd Jy(N_, N_);

        ode_system_.jacobian(t, y_base, dz_dt_eig, Jy, dZ_dt_sens_param);

        dZ_dt_sens.leftCols(N_).setZero();
        dZ_dt_sens += Jy * Z_sens;
      }

      inline void
      rhs_sens(const std::vector<stan::math::var>& initial,
               const std::vector<double>& param,
               const std::vector<double>& z,
               std::vector<double>& dz_dt,
               double t) const {
        using Eigen::Map;
        using Eigen::MatrixXd;
        using Eigen::VectorXd;

        const std::vector<double> y_base(z.begin(), z.begin() + N_);
        Map<VectorXd> dz_dt_eig(&dz_dt[0], N_);
        Map<MatrixXd> dZ_dt_sens(&dz_dt[0] + N_, N_, S_);
        Map<const MatrixXd> Z_sens(&z[0] + N_, N_, S_);

        MatrixXd Jy(N_, N_);

        ode_system_.jacobian(t, y_base, dz_dt_eig, Jy);

        dZ_dt_sens = Jy * Z_sens;
      }

      inline void
      rhs_sens(const std::vector<double>& initial,
               const std::vector<stan::math::var>& param,
               const std::vector<double>& z,
               std::vector<double>& dz_dt,
               double t) const {
        using Eigen::Map;
        using Eigen::MatrixXd;
        using Eigen::VectorXd;

        const std::vector<double> y_base(z.begin(), z.begin() + N_);
        Map<VectorXd> dz_dt_eig(&dz_dt[0], N_);
        Map<MatrixXd> dZ_dt_sens(&dz_dt[0] + N_, N_, S_);
        Map<const MatrixXd> Z_sens(&z[0] + N_, N_, S_);

        MatrixXd Jy(N_, N_);

        ode_system_.jacobian(t, y_base, dz_dt_eig, Jy, dZ_dt_sens);

        dZ_dt_sens += Jy * Z_sens;
      }
    };
  }  // math
}  // stan

#endif
