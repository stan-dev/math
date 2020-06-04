#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_STD_VECTOR_INTERFACE_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_STD_VECTOR_INTERFACE_ADAPTER_HPP

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
template <typename F>
struct integrate_ode_std_vector_interface_adapter {
  const F f_;
  const int num_vars__;

  integrate_ode_std_vector_interface_adapter(const F& f) : f_(f), num_vars__(f.num_vars__) {}

  template <typename T0, typename T1, typename T2>
  auto operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                  std::ostream* msgs, const std::vector<T2>& theta,
                  const std::vector<double>& x,
                  const std::vector<int>& x_int) const {
    return to_vector(f_(t, to_array_1d(y), msgs, theta, x, x_int));
  }
  template<typename v>
  void save_varis(v** p) const { f_.save_varis(p); }
  void accumulate_adjoints(double* p) const { f_.accumulate_adjoints(p); }
  void set_zer_adjoints() const { f_.set_zero_adjoints(); }

  struct DeepCopy_cl__ {
    const typename F::DeepCopy__ f_;
    const int num_vars__;
    struct ValueOf_cl__ {
      const typename F::ValueOf__ f_;
      const static int num_vars__ = 0;
      ValueOf_cl__(const integrate_ode_std_vector_interface_adapter<F>& f) : f_(f.f_) {}
      ValueOf_cl__(const DeepCopy_cl__& f) : f_(f.f_) {}
      template <typename T0, typename T1, typename T2>
      auto operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                      std::ostream* msgs, const std::vector<T2>& theta,
                      const std::vector<double>& x,
                      const std::vector<int>& x_int) const {
        return to_vector(f_(t, to_array_1d(y), msgs, theta, x, x_int));
      }
      template<typename v>
      void save_varis(v** p) const { f_.save_varis(p); }
      void accumulate_adjoints(double* p) const { f_.accumulate_adjoints(p); }
      void set_zer_adjoints() const { f_.set_zero_adjoints(); }
      using captured_scalar_t__ = double;
      using ValueOf__ = ValueOf_cl__;
      using DeepCopy__ = ValueOf_cl__;
    };
    DeepCopy_cl__(const integrate_ode_std_vector_interface_adapter<F>& f) : f_(f.f_), num_vars__(f.num_vars__) {}
    template <typename T0, typename T1, typename T2>
    auto operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                    std::ostream* msgs, const std::vector<T2>& theta,
                    const std::vector<double>& x,
                    const std::vector<int>& x_int) const {
      return to_vector(f_(t, to_array_1d(y), msgs, theta, x, x_int));
    }
    template<typename v>
    void save_varis(v** p) const { f_.save_varis(p); }
    void accumulate_adjoints(double* p) const { f_.accumulate_adjoints(p); }
    void set_zer_adjoints() const { f_.set_zero_adjoints(); }
    using captured_scalar_t__ = typename F::captured_scalar_t__;
    using ValueOf__ = ValueOf_cl__;
    using DeepCopy__ = DeepCopy_cl__;
  };
  using captured_scalar_t__ = typename F::captured_scalar_t__;
  using DeepCopy__ = DeepCopy_cl__;
  using ValueOf__ = typename DeepCopy__::ValueOf__;
};

}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
