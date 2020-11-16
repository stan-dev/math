#ifndef STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP

#include <stan/math/prim/meta/return_type.hpp>
#include <ostream>

namespace stan {
namespace math {

template <typename F>
struct closure_adapter {
  using captured_scalar_t__ = double;
  using ValueOf__ = closure_adapter<F>;
  static const size_t vars_count__ = 0;
  F f_;

  explicit closure_adapter(const F& f) : f_(f) {}

  template <typename... Args>
  auto operator()(std::ostream* msgs, Args... args) const {
    return f_(args..., msgs);
  }
  auto value_of__() const { return closure_adapter<F>(f_); }
  auto deep_copy_vars__() const { return closure_adapter<F>(f_); }
  void zero_adjoints__() const {}
  double* accumulate_adjoints__(double* dest) const { return dest; }
  template <typename Vari>
  Vari** save_varis(Vari** dest) const {
    return dest;
  }
};

template <typename F, typename T>
struct simple_closure {
  using captured_scalar_t__ = return_type_t<T>;
  using ValueOf__ = simple_closure<F, decltype(value_of(std::declval<T>()))>;
  const size_t vars_count__;
  F f_;
  T s_;

  explicit simple_closure(const F& f, T s)
      : f_(f), s_(s), vars_count__(count_vars(s)) {}

  template <typename... Args>
  auto operator()(std::ostream* msgs, Args... args) const {
    return f_(s_, args..., msgs);
  }
  auto value_of__() const { return ValueOf__(f_, value_of(s_)); }
  auto deep_copy_vars__() const {
    return simple_closure<F, T>(f_, deep_copy_vars(s_));
  }
  void zero_adjoints__() { zero_adjoints(s_); }
  double* accumulate_adjoints__(double* dest) const {
    return accumulate_adjoints(dest, s_);
  }
  template <typename Vari>
  Vari** save_varis__(Vari** dest) const {
    return save_varis(dest, s_);
  }
};

template <typename F>
auto from_lambda(F f) {
  return closure_adapter<F>(f);
}

template <typename F, typename T>
auto from_lambda(F f, T a) {
  return simple_closure<F, T>(f, a);
}

namespace internal {

template <typename F>
struct ode_closure_adapter {
  using captured_scalar_t__ = double;
  using ValueOf__ = ode_closure_adapter<F>;
  static const size_t vars_count__ = 0;
  const F f_;

  explicit ode_closure_adapter(const F& f) : f_(f) {}

  template <typename T0, typename T1, typename... Args>
  auto operator()(std::ostream* msgs, const T0& t,
                  const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                  Args... args) const {
    return f_(t, y, msgs, args...);
  }
  auto value_of__() const { return ode_closure_adapter<F>(f_); }
  auto deep_copy_vars__() const { return ode_closure_adapter<F>(f_); }
  void zero_adjoints__() const {}
  double* accumulate_adjoints__(double* dest) const { return dest; }
  template <typename Vari>
  Vari** save_varis(Vari** dest) const {
    return dest;
  }
};

template <typename F>
struct reduce_sum_closure_adapter {
  using captured_scalar_t__ = double;
  using ValueOf__ = reduce_sum_closure_adapter<F>;
  static const size_t vars_count__ = 0;
  const F f_;

  explicit reduce_sum_closure_adapter(const F& f) : f_(f) {}

  template <typename T, typename... Args>
  auto operator()(std::ostream* msgs, const std::vector<T>& sub_slice,
		  std::size_t start, std::size_t end,
                  Args... args) const {
    return f_(sub_slice, start, end, msgs, args...);
  }
  auto value_of__() const { return reduce_sum_closure_adapter<F>(f_); }
  auto deep_copy_vars__() const { return reduce_sum_closure_adapter<F>(f_); }
  void zero_adjoints__() const {}
  double* accumulate_adjoints__(double* dest) const { return dest; }
  template <typename Vari>
  Vari** save_varis(Vari** dest) const {
    return dest;
  }
};

}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
