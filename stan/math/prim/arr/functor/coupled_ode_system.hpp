#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_HPP

#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * The <code>coupled_ode_system</code> represents the coupled ode
 * system, which is the base ode and the sensitivities of the base ode
 * (derivatives with respect to the parameters of the base ode).
 *
 * This class provides a functor to be used by ode solvers, a class method
 * for the size of the coupled ode system, a class method to retrieve the
 * initial state, and a class method to convert from the coupled ode system
 * back to the base ode.
 *
 * @tparam F base ode system functor. Must provide
 *   <code>operator()(double t, std::vector<T1> y, std::vector<T2> theta,
 *          std::vector<double> x, std::vector<int>x_int, std::ostream*
 * msgs)</code>
 * @tparam T1 scalar type of the initial state
 * @tparam T2 scalar type of the parameters
 */
template <typename F, typename T1, typename T2>
struct coupled_ode_system {};

/**
 * <code>coupled_ode_system</code> specialization for for known
 * initial values and known parameters.
 *
 * For this case, the coupled ode is the same as the base ode. There
 * are no sensitivity parameters and the size of the coupled ode
 * system is the size of the base ode system.
 *
 * @tparam F base ode system functor. Must provide
 *   <code>operator()(double t, std::vector<double> y,
 *   std::vector<double> theta, std::vector<double> x,
 *   std::vector<int>x_int, std::ostream* msgs)</code>
 */
template <typename F>
class coupled_ode_system<F, double, double> {
 public:
  const F& f_;
  const std::vector<double>& y0_dbl_;
  const std::vector<double>& theta_dbl_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  const size_t N_;
  const size_t M_;
  const size_t size_;
  std::ostream* msgs_;

  /**
   * Construct the coupled ode system from the base system function,
   * initial state of the base system, parameters, data and a stream
   * for messages.
   *
   * @param[in] f base ode system functor
   * @param[in] y0 initial state of the base ode
   * @param[in] theta parameters of the base ode
   * @param[in] x real data
   * @param[in] x_int integer data
   * @param[in, out] msgs stream for messages
   */
  coupled_ode_system(const F& f, const std::vector<double>& y0,
                     const std::vector<double>& theta,
                     const std::vector<double>& x,
                     const std::vector<int>& x_int, std::ostream* msgs)
      : f_(f),
        y0_dbl_(y0),
        theta_dbl_(theta),
        x_(x),
        x_int_(x_int),
        N_(y0.size()),
        M_(theta.size()),
        size_(N_),
        msgs_(msgs) {}

  /**
   * Calculates the derivative of the coupled ode system with respect to time.
   *
   * @param[in] y current state of the coupled ode. This must be the
   *   correct size, <code>N_</code>.
   * @param[out] dy_dt populated with derivatives of the coupled
   *   system evaluated at specified state and time. This vector will be
   * overwritten.
   * @param[in] t time.
   * @throw exception if the system function does not return
   * a derivative vector of the same size as the state vector.
   */
  void operator()(const std::vector<double>& y, std::vector<double>& dy_dt,
                  double t) const {
    dy_dt = f_(t, y, theta_dbl_, x_, x_int_, msgs_);
    check_size_match("coupled_ode_system", "y", y.size(), "dy_dt",
                     dy_dt.size());
  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  int size() const { return size_; }

  /**
   * Returns the initial state of the coupled system. Here, it is
   * identical to base ode state because the initial state is known.
   *
   * @return initial state of the coupled system
   */
  std::vector<double> initial_state() const {
    std::vector<double> state(size_, 0.0);
    for (size_t n = 0; n < N_; n++) {
      state[n] = y0_dbl_[n];
    }
    return state;
  }
};

}  // namespace math
}  // namespace stan
#endif
