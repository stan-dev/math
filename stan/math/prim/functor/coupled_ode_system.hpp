#ifndef STAN_MATH_PRIM_FUNCTOR_COUPLED_ODE_SYSTEM_HPP
#define STAN_MATH_PRIM_FUNCTOR_COUPLED_ODE_SYSTEM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

template <bool ArithmeticArguments, typename F, typename T_y0, typename... Args>
struct coupled_ode_system_impl;

/**
 * The <code>coupled_ode_system_impl</code> for arithmetic arguments reduces to
 * the regular ode system (there are no sensitivities)
 *
 * @tparam F base ode system functor. Must provide
 *   <code>
 *     template<typename T_y, typename... T_args>
 *     operator()(double t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
 *     std::ostream* msgs, const T_args&... args)</code>
 */
template <typename F, typename T_y0, typename... Args>
struct coupled_ode_system_impl<true, F, T_y0, Args...> {
  const F& f_;
  const Eigen::VectorXd& y0_;
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
   * @param[in, out] msgs stream for messages
   * @param[in] args other additional arguments
   */
  coupled_ode_system_impl(const F& f, const Eigen::VectorXd& y0,
                          std::ostream* msgs, const Args&... args)
      : f_(f), y0_(y0), args_tuple_(args...), N_(y0.size()), msgs_(msgs) {}

  /**
   * Evaluate the right hand side of the coupled system at state z and time t,
   * and store the sensitivities in dz_dt.
   *
   * For all arithmetic types the coupled system right hand side is the same
   * as the regular ODE right hand side.
   *
   * @param[in] z State of coupled system
   * @param[out] dz_dt Evaluation of right hand side of coupled system
   * @param[in] t Time at which to evaluated right hand side
   * @throw exception if the system function does not return
   * the expected number of derivatives, N.
   */
  void operator()(const std::vector<double>& z, std::vector<double>& dz_dt,
                  double t) const {
    Eigen::VectorXd y(z.size());
    for (size_t n = 0; n < z.size(); ++n)
      y.coeffRef(n) = z[n];

    dz_dt.resize(y.size());

    Eigen::VectorXd f_y_t
        = apply([&](const Args&... args) { return f_(t, y, msgs_, args...); },
                args_tuple_);

    check_size_match("coupled_ode_system", "dy_dt", f_y_t.size(), "states",
                     y.size());

    Eigen::Map<Eigen::VectorXd>(dz_dt.data(), dz_dt.size()) = f_y_t;
  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  size_t size() const { return N_; }

  /**
   * Returns the initial state of the coupled system. For arithmetic types
   * this is the same as the regular ode system.
   *
   * @return the initial condition of the coupled system
   */
  std::vector<double> initial_state() const {
    std::vector<double> initial(size(), 0.0);

    for (size_t i = 0; i < N_; i++) {
      initial[i] = value_of(y0_(i));
    }

    return initial;
  }
};

template <typename F, typename T_y0, typename... Args>
struct coupled_ode_system
    : public coupled_ode_system_impl<
          std::is_arithmetic<return_type_t<T_y0, Args...>>::value, F, T_y0,
          Args...> {
  coupled_ode_system(const F& f,
                     const Eigen::Matrix<T_y0, Eigen::Dynamic, 1>& y0,
                     std::ostream* msgs, const Args&... args)
      : coupled_ode_system_impl<
          std::is_arithmetic<return_type_t<T_y0, Args...>>::value, F, T_y0,
          Args...>(f, y0, msgs, args...) {}
};

}  // namespace math
}  // namespace stan
#endif
