#ifndef STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP

#include <stan/math/prim/meta/error_index.hpp>
#include <stan/math/prim/meta/is_stan_closure.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <ostream>

namespace stan {
namespace math {
namespace internal {

/**
 * A closure that wraps a C++ lambda and captures values.
 */
template <bool Ref, typename F, typename... Ts>
struct base_closure {
  using captured_scalar_t__ = return_type_t<Ts...>;
  using ValueOf__
      = base_closure<false, F,
                         decltype(eval(value_of(std::declval<Ts>())))...>;
  using CopyOf__ = base_closure<false, F, Ts...>;
  F f_;
  std::tuple<capture_type_t<Ts, Ref>...> captures_;

  explicit base_closure(const F& f, const Ts&... args)
      : f_(f), captures_(args...) {}

  template <typename... Args>
  auto operator()(std::ostream* msgs, Args... args) const {
    return apply(
        [this, msgs, args...](auto... s) { return f_(s..., args..., msgs); },
        captures_);
  }
  size_t count_vars__() const {
    return apply([this](auto... s) { return count_vars(s...); }, captures_);
  }
  auto value_of__() const {
    return apply(
        [this](auto... s) { return ValueOf__(f_, eval(value_of(s))...); },
        captures_);
  }
  auto deep_copy_vars__() const {
    return apply(
        [this](auto... s) { return CopyOf__(f_, eval(deep_copy_vars(s))...); },
        captures_);
  }
  void zero_adjoints__() {
    apply([](auto... s) { zero_adjoints(s...); }, captures_);
  }
  double* accumulate_adjoints__(double* dest) const {
    return apply([dest](auto... s) { return accumulate_adjoints(dest, s...); },
                 captures_);
  }
  template <typename Vari>
  Vari** save_varis__(Vari** dest) const {
    return apply([dest](auto... s) { return save_varis(dest, s...); },
                 captures_);
  }
};

/**
 * A closure that takes rng argument.
 */
template <typename F>
struct empty_closure_rng {
  using captured_scalar_t__ = double;
  using ValueOf__ = empty_closure_rng<F>;
  using CopyOf__ = empty_closure_rng<F>;
  F f_;

  explicit empty_closure_rng(const F& f) : f_(f) {}

  template <typename Rng, typename... Args>
  auto operator()(const Rng& rng, std::ostream* msgs, Args... args) const {
    return f_(args..., rng, msgs);
  }
  size_t count_vars__() const { return 0; }
  auto value_of__() const { return ValueOf__(f_); }
  auto copy_of__() const { return CopyOf__(f_); }
  auto deep_copy_vars__() const { return CopyOf__(f_); }
  void zero_adjoints__() const {}
  double* accumulate_adjoints__(double* dest) const { return dest; }
  template <typename Vari>
  Vari** save_varis(Vari** dest) const {
    return dest;
  }
};

/**
 * A closure that can be called with `propto` template argument.
 */
template <typename F>
struct empty_closure_lpdf {
  using captured_scalar_t__ = double;
  using ValueOf__ = empty_closure_lpdf<F>;
  using CopyOf__ = empty_closure_lpdf<F>;
  F f_;

  explicit empty_closure_lpdf(const F& f) : f_(f) {}

  template <bool propto = false, typename... Args>
  auto operator()(std::ostream* msgs, Args... args) const {
    return f_.template operator()<propto>(args..., msgs);
  }
  size_t count_vars__() const { return 0; }
  auto value_of__() const { return ValueOf__(f_); }
  auto copy_of__() const { return CopyOf__(f_); }
  auto deep_copy_vars__() const { return CopyOf__(f_); }
  void zero_adjoints__() const {}
  double* accumulate_adjoints__(double* dest) const { return dest; }
  template <typename Vari>
  Vari** save_varis(Vari** dest) const {
    return dest;
  }
};

/**
 * A closure that accesses logprob accumulator.
 */
template <typename F>
struct empty_closure_lp {
  using captured_scalar_t__ = double;
  using ValueOf__ = empty_closure_lp<F>;
  using CopyOf__ = empty_closure_lp<F>;
  static const size_t vars_count__ = 0;
  F f_;

  explicit empty_closure_lp(const F& f) : f_(f) {}

  template <typename T_lp, typename T_lp_accum, typename... Args>
  auto operator()(T_lp_accum& lp, T_lp& lp_accum, std::ostream* msgs,
                  Args... args) const {
    return f_(args..., lp, lp_accum, msgs);
  }
  size_t count_vars__() const { return 0; }
  auto value_of__() const { return ValueOf__(f_); }
  auto copy_of__() const { return CopyOf__(f_); }
  auto deep_copy_vars__() const { return CopyOf__(f_); }
  void zero_adjoints__() const {}
  double* accumulate_adjoints__(double* dest) const { return dest; }
  template <typename Vari>
  Vari** save_varis(Vari** dest) const {
    return dest;
  }
};

}  // namespace internal

/**
 * Higher-order functor suitable for calling a closure inside variadic ODE
 * solvers.
 */
struct ode_closure_adapter {
  template <typename F, typename T0, typename T1, typename... Args>
  auto operator()(const T0& t, const T1& y, std::ostream* msgs, const F& f,
                  Args... args) const {
    return f(msgs, t, y, args...);
  }
};

/**
 * Create a closure from a C++ lambda and captures.
 */
template <typename F, typename... Ts>
auto from_lambda(const F& f, const Ts&... a) {
  return internal::base_closure<true, F, Ts...>(f, a...);
}

/**
 * Create a closure from an rng functor.
 */
template <typename F>
auto rng_from_lambda(const F& f) {
  return internal::empty_closure_rng<F>(f);
}

/**
 * Create a closure from an lpdf functor.
 */
template <typename F>
auto lpdf_from_lambda(const F& f) {
  return internal::empty_closure_lpdf<F>(f);
}

/**
 * Create a closure from a functor that needs access to logprob accumulator.
 */
template <typename F>
auto lp_from_lambda(const F& f) {
  return internal::empty_closure_lp<F>(f);
}

/**
 * A wrapper that sets propto template argument when calling the inner closure.
 */
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

/**
 * Higher-order functor that invokes a closure inside a reduce_sum call.
 */
struct reduce_sum_closure_adapter {
  template <typename F, typename T, typename... Args>
  auto operator()(const std::vector<T>& sub_slice, std::size_t start,
                  std::size_t end, std::ostream* msgs, const F& f,
                  Args... args) const {
    return f(msgs, sub_slice, start + error_index::value,
             end + error_index::value, args...);
  }
};

}  // namespace math
}  // namespace stan

#endif
