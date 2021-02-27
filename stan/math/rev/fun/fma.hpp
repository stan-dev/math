#ifndef STAN_MATH_REV_FUN_FMA_HPP
#define STAN_MATH_REV_FUN_FMA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fma.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

/**
 * The fused multiply-add function for three variables (C99).
 * This function returns the product of the first two arguments
 * plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
inline var fma(const var& a, const var& b, const var& c) {
  return make_callback_var(fma(a.val(), b.val(), c.val()), [a, b, c](auto& vi) {
    a.adj() += vi.adj() * b.val();
    b.adj() += vi.adj() * a.val();
    c.adj() += vi.adj();
  });
}

/**
 * The fused multiply-add function for two variables and a value
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + c = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + c = x\f$.
 *
 * @tparam Tc type of the summand
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tc, require_arithmetic_t<Tc>* = nullptr>
inline var fma(const var& a, const var& b, Tc&& c) {
  return make_callback_var(fma(a.val(), b.val(), c), [a, b](auto& vi) {
    a.adj() += vi.adj() * b.val();
    b.adj() += vi.adj() * a.val();
  });
}

/**
 * The fused multiply-add function for a variable, value, and
 * variable (C99).  This function returns the product of the first
 * two arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * c) + z = c\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * c) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tb type of the second multiplicand
 * @tparam Tc type of the summand
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tb, typename Tc,
          require_arithmetic_t<Tb>* = nullptr,
          require_all_var_t<Ta, Tc>* = nullptr>
inline var fma(Ta&& a, Tb&& b, Tc&& c) {
  return make_callback_var(fma(a.val(), b, c.val()), [a, b, c](auto& vi) {
    a.adj() += vi.adj() * b;
    c.adj() += vi.adj();
  });
}

/**
 * The fused multiply-add function for a variable and two values
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The double-based version
 * <code>::%fma(double, double, double)</code> is defined in
 * <code>&lt;cmath&gt;</code>.
 *
 * The derivative is
 *
 * \f$\frac{d}{d x} (x * c) + d = c\f$.
 *
 * @tparam Tb type of the second multiplicand
 * @tparam Tc type of the summand
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tb, typename Tc, require_all_arithmetic_t<Tb, Tc>* = nullptr>
inline var fma(const var& a, Tb&& b, Tc&& c) {
  return make_callback_var(fma(a.val(), b, c), [a, b](auto& vi) {
    a.adj() += vi.adj() * b;
  });
}

/**
 * The fused multiply-add function for a value, variable, and
 * value (C99).  This function returns the product of the first
 * two arguments plus the third argument.
 *
 * The derivative is
 *
 * \f$\frac{d}{d y} (c * y) + d = c\f$, and
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tc type of the summand
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tc, require_all_arithmetic_t<Ta, Tc>* = nullptr>
inline var fma(Ta&& a, const var& b, Tc&& c) {
  return make_callback_var(fma(a, b.val(), c), [a, b](auto& vi) {
    b.adj() += vi.adj() * a;
  });
}

/**
 * The fused multiply-add function for two values and a variable,
 * and value (C99).  This function returns the product of the
 * first two arguments plus the third argument.
 *
 * The derivative is
 *
 * \f$\frac{\partial}{\partial z} (c * d) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tb type of the second multiplicand
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tb, require_all_arithmetic_t<Ta, Tb>* = nullptr>
inline var fma(Ta&& a, Tb&& b, const var& c) {
  return make_callback_var(fma(a, b, c.val()), [c](auto& vi) {
    c.adj() += vi.adj();
  });
}

/**
 * The fused multiply-add function for a value and two variables
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial y} (c * y) + z = c\f$, and
 *
 * \f$\frac{\partial}{\partial z} (c * y) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, require_arithmetic_t<Ta>* = nullptr>
inline var fma(Ta&& a, const var& b, const var& c) {
  return make_callback_var(fma(a, b.val(), c.val()), [a, b, c](auto& vi) {
    b.adj() += vi.adj() * a;
    c.adj() += vi.adj();
  });
}

namespace internal {
/**
 * Overload for matrix, matrix, matrix
 */
template <typename T1, typename T2, typename T3, typename T4,
 require_all_matrix_t<T1, T2, T3>* = nullptr>
inline auto fma_reverse_pass(T1& arena_a, T2& arena_b, T3& arena_c, T4& ret) {
  return [arena_a, arena_b, arena_c, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_a).adj().array() += ret.adj().array() * value_of(arena_b).array();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_b).adj().array() += ret.adj().array() * value_of(arena_a).array();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_c).adj().array() += ret.adj().array();
    }
  };
}

/**
 * Overload for scalar, matrix, matrix
 */
template <typename T1, typename T2, typename T3, typename T4,
 require_all_matrix_t<T2, T3>* = nullptr,
 require_stan_scalar_t<T1>* = nullptr>
inline auto fma_reverse_pass(T1& arena_a, T2& arena_b, T3& arena_c, T4& ret) {
  return [arena_a, arena_b, arena_c, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_a).adj() += (ret.adj().array() * value_of(arena_b).array()).sum();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_b).adj().array() += ret.adj().array() * value_of(arena_a);
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_c).adj().array() += ret.adj().array();
    }
  };
}

/**
 * Overload for matrix, scalar, matrix
 */
 template <typename T1, typename T2, typename T3, typename T4,
 require_all_matrix_t<T1, T3>* = nullptr,
 require_stan_scalar_t<T2>* = nullptr>
inline auto fma_reverse_pass(T1& arena_a, T2& arena_b, T3& arena_c, T4& ret) {
  return [arena_a, arena_b, arena_c, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_a).adj().array() += ret.adj().array() * value_of(arena_b);
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_b).adj() += (ret.adj().array() * value_of(arena_a).array()).sum();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_c).adj().array() += ret.adj().array();
    }
  };
}

/**
 * Overload for scalar, scalar, matrix
 */
 template <typename T1, typename T2, typename T3, typename T4,
 require_matrix_t<T3>* = nullptr,
 require_all_stan_scalar_t<T1, T2>* = nullptr>
inline auto fma_reverse_pass(T1& arena_a, T2& arena_b, T3& arena_c, T4& ret) {
  return [arena_a, arena_b, arena_c, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_a).adj() += (ret.adj().array() * value_of(arena_b)).sum();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_b).adj() += (ret.adj().array() * value_of(arena_a)).sum();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_c).adj().array() += ret.adj().array();
    }
  };
}


/**
 * Overload for matrix, matrix, scalar
 */
 template <typename T1, typename T2, typename T3, typename T4,
 require_all_matrix_t<T1, T2>* = nullptr,
 require_stan_scalar_t<T3>* = nullptr>
inline auto fma_reverse_pass(T1& arena_a, T2& arena_b, T3& arena_c, T4& ret) {
  return [arena_a, arena_b, arena_c, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_a).adj().array() += ret.adj().array() * value_of(arena_b).array();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_b).adj().array() += ret.adj().array() * value_of(arena_a).array();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_c).adj() += ret.adj().sum();
    }
  };
}

/**
 * Overload for scalar, matrix, scalar
 */
 template <typename T1, typename T2, typename T3, typename T4,
 require_matrix_t<T2>* = nullptr,
 require_all_stan_scalar_t<T1, T3>* = nullptr>
inline auto fma_reverse_pass(T1& arena_a, T2& arena_b, T3& arena_c, T4& ret) {
  return [arena_a, arena_b, arena_c, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_a).adj() += (ret.adj().array() * value_of(arena_b).array()).sum();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_b).adj().array() += ret.adj().array() * value_of(arena_a);
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_c).adj() += ret.adj().sum();
    }
  };
}


/**
 * Overload for matrix, scalar, scalar
 */
template <typename T1, typename T2, typename T3, typename T4,
 require_matrix_t<T1>* = nullptr,
 require_all_stan_scalar_t<T2, T3>* = nullptr>
inline auto fma_reverse_pass(T1& arena_a, T2& arena_b, T3& arena_c, T4& ret) {
  return [arena_a, arena_b, arena_c, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_a).adj().array() += ret.adj().array() * value_of(arena_b);
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_b).adj() += (ret.adj().array() * value_of(arena_a).array()).sum();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_c).adj() += ret.adj().sum();
    }
  };
}

}

/**
 * The fused multiply-add function for three variables (C99).
 * This function returns the product of the first two arguments
 * plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename T1, typename T2, typename T3,
 require_any_matrix_t<T1, T2, T3>* = nullptr,
 require_var_t<return_type_t<T1, T2, T3>>* = nullptr>
inline auto fma(const T1& a, const T2& b, const T3& c) {
  arena_t<T1> arena_a = a;
  arena_t<T2> arena_b = b;
  arena_t<T3> arena_c = c;
  using inner_ret_type = decltype(fma(value_of(arena_a), value_of(arena_b), value_of(arena_c)));
  using ret_type = return_var_matrix_t<inner_ret_type, T1, T2, T3>;
  arena_t<ret_type> ret = fma(value_of(arena_a), value_of(arena_b), value_of(arena_c));
  reverse_pass_callback(internal::fma_reverse_pass(arena_a, arena_b, arena_c, ret));
  return ret_type(ret);
}


}  // namespace math
}  // namespace stan
#endif
