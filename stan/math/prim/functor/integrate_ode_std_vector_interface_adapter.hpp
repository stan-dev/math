#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_STD_VECTOR_INTERFACE_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_STD_VECTOR_INTERFACE_ADAPTER_HPP

#include <stan/math/prim/fun/to_array_1d.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

namespace internal {

template <typename F>
struct integrate_ode_std_vector_interface_adapter {
  // NOTE: I made this copy by value because
  //  this function might go out of scope in reverse mode
  const F f_;

  integrate_ode_std_vector_interface_adapter(const F& f) : f_(f) {
  }
  
  template<typename T0, typename T1, typename T2>
  auto operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y, std::ostream* msgs,
		  const std::vector<T2>& theta, const std::vector<double>& x,
		  const std::vector<int>& x_int) const {
    return to_vector(f_(t, to_array_1d(y), msgs, theta, x, x_int));
  }
};

}

}  // namespace math
}  // namespace stan

#endif
