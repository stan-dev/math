#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_1D_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_1D_ADAPTER_HPP

#include <ostream>
#include <vector>

/**
 * Adapt the non-variadic integrate_1d arguments to the variadic
 * integrate_1d_impl interface
 *
 * @tparam F type of function to adapt
 */
template <typename F>
struct integrate_1d_adapter {
  const F& f_;

  explicit integrate_1d_adapter(const F& f) : f_(f) {}

  template <typename T_a, typename T_b, typename T_theta>
  auto operator()(const T_a& x, const T_b& xc, std::ostream* msgs,
                  const std::vector<T_theta>& theta,
                  const std::vector<double>& x_r,
                  const std::vector<int>& x_i) const {
    return f_(x, xc, theta, x_r, x_i, msgs);
  }
};

#endif
