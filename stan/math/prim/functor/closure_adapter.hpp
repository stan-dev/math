#ifndef STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP

#include <stan/math/prim/meta/return_type.hpp>
#include <ostream>

namespace stan {
namespace math {

template <typename F>
struct empty_closure {
  using captured_scalar_t__ = double;
  using ValueOf__ = empty_closure<F>;
  using CopyOf__ = empty_closure<F>;
  static const size_t vars_count__ = 0;
  F f_;

  explicit empty_closure(const F& f) : f_(f) {}

  template <typename... Args>
  auto operator()(std::ostream* msgs, Args... args) const {
    return f_(args..., msgs);
  }
  size_t count_vars__() const { return 0; }
  auto value_of__() const { return empty_closure<F>(f_); }
  auto shallow_copy__() const { return empty_closure<F>(f_); }
  auto deep_copy_vars__() const { return empty_closure<F>(f_); }
  void zero_adjoints__() const {}
  double* accumulate_adjoints__(double* dest) const { return dest; }
  template <typename Vari>
  Vari** save_varis(Vari** dest) const {
    return dest;
  }
};

template <typename F, typename T>
struct one_arg_closure {
  using captured_scalar_t__ = return_type_t<T>;
  using ValueOf__ = one_arg_closure<F, decltype(value_of(std::declval<T>()))>;
  using CopyOf__ = one_arg_closure<F, T>;
  F f_;
  T s_;

  explicit one_arg_closure(const F& f, const T& s) : f_(f), s_(s) {}

  template <typename... Args>
  auto operator()(std::ostream* msgs, Args... args) const {
    return f_(s_, args..., msgs);
  }
  size_t count_vars__() const { return count_vars(s_); }
  auto value_of__() const { return ValueOf__(f_, value_of(s_)); }
  auto shallow_copy__() const { return one_arg_closure<F, T>(f_, s_); }
  auto deep_copy_vars__() const {
    return one_arg_closure<F, T>(f_, deep_copy_vars(s_));
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
  return empty_closure<F>(f);
}

template <typename F, typename T>
auto from_lambda(F f, T a) {
  return one_arg_closure<F, T>(f, a);
}

template <bool Propto, typename F, bool Ref>
struct lpdf_wrapper {
  using captured_scalar_t__ = return_type_t<F>;
  using ValueOf__
      = lpdf_wrapper<Propto, decltype(std::declval<F>().value_of__()), false>;
  using CopyOf__
      = lpdf_wrapper<Propto, decltype(std::declval<F>().copy_of__()), false>;
  capture_type_t<F, Ref> f_;

  explicit lpdf_wrapper(const F& f) : f_(f) {}

  template <bool propto>
  auto with_propto() {
    return lpdf_wrapper < Propto && propto, F, true > (f_);
  }

  template <bool propto = Propto, typename... Args>
  auto operator()(Args... args) const {
    return f_.template operator() < Propto && propto > (args...);
  }
  size_t count_vars__() const { return count_vars(f_); }
  auto value_of__() const { return ValueOf__(value_of(f_)); }
  auto deep_copy_vars__() const { return CopyOf__(deep_copy_vars(f_)); }
  auto copy_of__() const { return CopyOf__(f_.copy_of__()); }
  void zero_adjoints__() { zero_adjoints(f_); }
  double* accumulate_adjoints__(double* dest) const {
    return accumulate_adjoints(dest, f_);
  }
  template <typename Vari>
  Vari** save_varis__(Vari** dest) const {
    return save_varis(dest, f_);
  }
};

struct reduce_sum_closure_adapter {
  template <typename F, typename T, typename... Args>
  auto operator()(const std::vector<T>& sub_slice, std::size_t start,
                  std::size_t end, std::ostream* msgs, const F& f,
                  Args... args) const {
    return f(msgs, sub_slice, start, end, args...);
  }
};

namespace internal {

struct ode_closure_adapter {
  template <typename F, typename T0, typename T1, typename... Args>
  auto operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                  std::ostream* msgs, const F& f, Args... args) const {
    return f(msgs, t, y, args...);
  }
};

}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
