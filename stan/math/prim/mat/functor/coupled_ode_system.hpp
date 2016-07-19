#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_COUPLED_ODE_SYSTEM_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_COUPLED_ODE_SYSTEM_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
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
    class coupled_ode_system;

    template <typename F>
    class coupled_ode_system<F, double, double> {
      const F& f_;
      const std::vector<double>& y0_;
      const std::vector<double>& theta_;
      const size_t N_;
      const size_t M_;
      const size_t S_;
      const size_t size_;
      const std::vector<double>& x_;
      const std::vector<int>& x_int_;
      std::ostream* msgs_;

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
                         const std::vector<double>& y0,
                         const std::vector<double>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         std::ostream* msgs)
        : f_(f),
          y0_(y0),
          theta_(theta),
          N_(y0.size()),
          M_(theta.size()),
          S_(0),
          size_(N_ * (1 + S_)),
          x_(x),
          x_int_(x_int),
          msgs_(msgs) { }

      /**
       * Returns the ODE RHS.
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
        dz_dt = f_(t, z, theta_, x_, x_int_, msgs_);
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
       * Returns the initial state of the system.
       *
       * @return the initial condition of the base system.  This is
       * a vector of length size().
       */
      std::vector<double> initial_state() const {
        return y0_;
      }
    };
  }  // math
}  // stan

#endif
