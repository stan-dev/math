#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_STD_VECTOR_INTERFACE_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_STD_VECTOR_INTERFACE_ADAPTER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/to_array_1d.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

namespace internal {

/**
 * Wrap a functor designed for use with integrate_ode_bdf, integrate_ode_rk45,
 * and integrate_ode_adams to use with the new ode_bdf/ode_rk45 interfaces.
 *
 * The old functors took the ODE state as a std::vector. The new ones take
 * state as an Eigen::Matrix. The adapter converts to and from these forms
 * so that the old ODE interfaces can work.
 */
template <bool is_closure, typename F>
struct integrate_ode_std_vector_interface_adapter_impl;

template <typename F>
struct integrate_ode_std_vector_interface_adapter_impl<false, F> {
  const F& f_;
  explicit integrate_ode_std_vector_interface_adapter_impl(const F& f)
      : f_(f) {}

  template <typename T0, typename T1, typename T2>
  auto operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                  std::ostream* msgs, const std::vector<T2>& theta,
                  const std::vector<double>& x,
                  const std::vector<int>& x_int) const {
    return to_vector(f_(t, to_array_1d(y), theta, x, x_int, msgs));
  }
};

template <typename F>
struct integrate_ode_std_vector_interface_adapter_impl<true, F> {
  using captured_scalar_t__ = typename F::captured_scalar_t__;
  using ValueOf__
      = integrate_ode_std_vector_interface_adapter_impl<true,
                                                        typename F::ValueOf__>;
  F f_;

  explicit integrate_ode_std_vector_interface_adapter_impl(const F& f)
      : f_(f) {}

  template <typename T0, typename T1, typename T2>
  auto operator()(std::ostream* msgs, const T0& t,
                  const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                  const std::vector<T2>& theta, const std::vector<double>& x,
                  const std::vector<int>& x_int) const {
    return to_vector(f_(msgs, t, to_array_1d(y), theta, x, x_int));
  }

  size_t count_vars__() const { return f_.count_vars__(); }
  auto value_of__() const { return ValueOf__(f_.value_of__()); }
  auto deep_copy_vars__() const {
    return integrate_ode_std_vector_interface_adapter_impl<true, F>(
        f_.deep_copy_vars__());
  }
  void zero_adjoints__() { f_.zero_adjoints__(); }
  double* accumulate_adjoints__(double* dest) const {
    return f_.accumulate_adjoints__(dest);
  }
  template <typename Vari>
  Vari** save_varis__(Vari** dest) const {
    return f_.save_varis__(dest);
  }
};

template <typename F>
using integrate_ode_std_vector_interface_adapter
    = integrate_ode_std_vector_interface_adapter_impl<is_stan_closure<F>::value,
                                                      F>;

}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
