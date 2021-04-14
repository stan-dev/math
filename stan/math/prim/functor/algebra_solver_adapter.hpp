#ifndef STAN_MATH_PRIM_FUNCTOR_ALGEBRA_SOLVER_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_ALGEBRA_SOLVER_ADAPTER_HPP

#include <ostream>
#include <vector>

/**
 * Adapt the non-variadic algebra_solver_newton and algebra_solver_powell
 * arguemts to the variadic algebra_solver_newton_impl and
 * algebra_solver_powell_impl interfaces.
 *
 * @tparam F type of function to adapt.
 */
template <typename F>
struct algebra_solver_adapter {
  const F& f_;

  explicit algebra_solver_adapter(const F& f) : f_(f) {}

  template <typename T1, typename T2, typename T3, typename T4>
  auto operator()(const T1& x, std::ostream* msgs, const T2& y, const T3& dat,
                  const T4& dat_int) const {
    return f_(x, y, dat, dat_int, msgs);
  }
};

#endif
