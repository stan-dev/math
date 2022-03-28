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
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
inline var fma(const var& x, const var& y, const var& z) {
  return make_callback_var(fma(x.val(), y.val(), z.val()), [x, y, z](auto& vi) {
    x.adj() += vi.adj() * y.val();
    y.adj() += vi.adj() * x.val();
    z.adj() += vi.adj();
  });
}

/**
 * The fused multiply-add function for two variables and a value
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$.
 *
 * @tparam Tc type of the summand
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tc, require_arithmetic_t<Tc>* = nullptr>
inline var fma(const var& x, const var& y, Tc&& z) {
  return make_callback_var(fma(x.val(), y.val(), z), [x, y](auto& vi) {
    x.adj() += vi.adj() * y.val();
    y.adj() += vi.adj() * x.val();
  });
}

/**
 * The fused multiply-add function for a variable, value, and
 * variable (C99).  This function returns the product of the first
 * two arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tb type of the second multiplicand
 * @tparam Tc type of the summand
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tb, require_arithmetic_t<Tb>* = nullptr>
inline var fma(const var& x, Tb&& y, const var& z) {
  return make_callback_var(fma(x.val(), y, z.val()), [x, y, z](auto& vi) {
    x.adj() += vi.adj() * y;
    z.adj() += vi.adj();
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
 * \f$\frac{d}{d x} (x * y) + z = y\f$.
 *
 * @tparam Tb type of the second multiplicand
 * @tparam Tc type of the summand
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tb, typename Tc, require_all_arithmetic_t<Tb, Tc>* = nullptr>
inline var fma(const var& x, Tb&& y, Tc&& z) {
  return make_callback_var(fma(x.val(), y, z),
                           [x, y](auto& vi) { x.adj() += vi.adj() * y; });
}

/**
 * The fused multiply-add function for a value, variable, and
 * value (C99).  This function returns the product of the first
 * two arguments plus the third argument.
 *
 * The derivative is
 *
 * \f$\frac{d}{d y} (x * y) + z = x\f$, and
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tc type of the summand
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tc, require_all_arithmetic_t<Ta, Tc>* = nullptr>
inline var fma(Ta&& x, const var& y, Tc&& z) {
  return make_callback_var(fma(x, y.val(), z),
                           [x, y](auto& vi) { y.adj() += vi.adj() * x; });
}

/**
 * The fused multiply-add function for two values and a variable,
 * and value (C99).  This function returns the product of the
 * first two arguments plus the third argument.
 *
 * The derivative is
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tb type of the second multiplicand
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tb, require_all_arithmetic_t<Ta, Tb>* = nullptr>
inline var fma(Ta&& x, Tb&& y, const var& z) {
  return make_callback_var(fma(x, y, z.val()),
                           [z](auto& vi) { z.adj() += vi.adj(); });
}

/**
 * The fused multiply-add function for a value and two variables
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, require_arithmetic_t<Ta>* = nullptr>
inline var fma(Ta&& x, const var& y, const var& z) {
  return make_callback_var(fma(x, y.val(), z.val()), [x, y, z](auto& vi) {
    y.adj() += vi.adj() * x;
    z.adj() += vi.adj();
  });
}

namespace internal {
/**
 * Overload for matrix, matrix, matrix
 */
template <typename T1, typename T2, typename T3, typename T4,
          require_all_matrix_t<T1, T2, T3>* = nullptr>
inline auto fma_reverse_pass(T1& arena_x, T2& arena_y, T3& arena_z, T4& ret) {
  return [arena_x, arena_y, arena_z, ret]() mutable {
    using T1_var = arena_t<plain_type_t<promote_scalar_t<var, T1>>>;
    using T2_var = arena_t<plain_type_t<promote_scalar_t<var, T2>>>;
    using T3_var = arena_t<plain_type_t<promote_scalar_t<var, T3>>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_x).adj().array()
          += ret.adj().array() * value_of(arena_y).array();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_y).adj().array()
          += ret.adj().array() * value_of(arena_x).array();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_z).adj().array() += ret.adj().array();
    }
  };
}

/**
 * Overload for scalar, matrix, matrix
 */
template <typename T1, typename T2, typename T3, typename T4,
          require_all_matrix_t<T2, T3>* = nullptr,
          require_stan_scalar_t<T1>* = nullptr>
inline auto fma_reverse_pass(T1& arena_x, T2& arena_y, T3& arena_z, T4& ret) {
  return [arena_x, arena_y, arena_z, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_x).adj()
          += (ret.adj().array() * value_of(arena_y).array()).sum();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_y).adj().array()
          += ret.adj().array() * value_of(arena_x);
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_z).adj().array() += ret.adj().array();
    }
  };
}

/**
 * Overload for matrix, scalar, matrix
 */
template <typename T1, typename T2, typename T3, typename T4,
          require_all_matrix_t<T1, T3>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto fma_reverse_pass(T1& arena_x, T2& arena_y, T3& arena_z, T4& ret) {
  return [arena_x, arena_y, arena_z, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_x).adj().array()
          += ret.adj().array() * value_of(arena_y);
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_y).adj()
          += (ret.adj().array() * value_of(arena_x).array()).sum();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_z).adj().array() += ret.adj().array();
    }
  };
}

/**
 * Overload for scalar, scalar, matrix
 */
template <typename T1, typename T2, typename T3, typename T4,
          require_matrix_t<T3>* = nullptr,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline auto fma_reverse_pass(T1& arena_x, T2& arena_y, T3& arena_z, T4& ret) {
  return [arena_x, arena_y, arena_z, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_x).adj()
          += (ret.adj().array() * value_of(arena_y)).sum();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_y).adj()
          += (ret.adj().array() * value_of(arena_x)).sum();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_z).adj().array() += ret.adj().array();
    }
  };
}

/**
 * Overload for matrix, matrix, scalar
 */
template <typename T1, typename T2, typename T3, typename T4,
          require_all_matrix_t<T1, T2>* = nullptr,
          require_stan_scalar_t<T3>* = nullptr>
inline auto fma_reverse_pass(T1& arena_x, T2& arena_y, T3& arena_z, T4& ret) {
  return [arena_x, arena_y, arena_z, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_x).adj().array()
          += ret.adj().array() * value_of(arena_y).array();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_y).adj().array()
          += ret.adj().array() * value_of(arena_x).array();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_z).adj() += ret.adj().sum();
    }
  };
}

/**
 * Overload for scalar, matrix, scalar
 */
template <typename T1, typename T2, typename T3, typename T4,
          require_matrix_t<T2>* = nullptr,
          require_all_stan_scalar_t<T1, T3>* = nullptr>
inline auto fma_reverse_pass(T1& arena_x, T2& arena_y, T3& arena_z, T4& ret) {
  return [arena_x, arena_y, arena_z, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_x).adj()
          += (ret.adj().array() * value_of(arena_y).array()).sum();
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_y).adj().array()
          += ret.adj().array() * value_of(arena_x);
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_z).adj() += ret.adj().sum();
    }
  };
}

/**
 * Overload for matrix, scalar, scalar
 */
template <typename T1, typename T2, typename T3, typename T4,
          require_matrix_t<T1>* = nullptr,
          require_all_stan_scalar_t<T2, T3>* = nullptr>
inline auto fma_reverse_pass(T1& arena_x, T2& arena_y, T3& arena_z, T4& ret) {
  return [arena_x, arena_y, arena_z, ret]() mutable {
    using T1_var = arena_t<promote_scalar_t<var, T1>>;
    using T2_var = arena_t<promote_scalar_t<var, T2>>;
    using T3_var = arena_t<promote_scalar_t<var, T3>>;
    if (!is_constant<T1>::value) {
      forward_as<T1_var>(arena_x).adj().array()
          += ret.adj().array() * value_of(arena_y);
    }
    if (!is_constant<T2>::value) {
      forward_as<T2_var>(arena_y).adj()
          += (ret.adj().array() * value_of(arena_x).array()).sum();
    }
    if (!is_constant<T3>::value) {
      forward_as<T3_var>(arena_z).adj() += ret.adj().sum();
    }
  };
}

}  // namespace internal

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
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename T1, typename T2, typename T3,
          require_any_matrix_t<T1, T2, T3>* = nullptr,
          require_var_t<return_type_t<T1, T2, T3>>* = nullptr>
inline auto fma(const T1& x, const T2& y, const T3& z) {
  arena_t<T1> arena_x = x;
  arena_t<T2> arena_y = y;
  arena_t<T3> arena_z = z;
  if (is_matrix<T1>::value && is_matrix<T2>::value) {
    check_matching_dims("fma", "x", arena_x, "y", arena_y);
  }
  if (is_matrix<T1>::value && is_matrix<T3>::value) {
    check_matching_dims("fma", "x", arena_x, "z", arena_z);
  }
  if (is_matrix<T2>::value && is_matrix<T3>::value) {
    check_matching_dims("fma", "y", arena_y, "z", arena_z);
  }
  using inner_ret_type
      = decltype(fma(value_of(arena_x), value_of(arena_y), value_of(arena_z)));
  using ret_type = return_var_matrix_t<inner_ret_type, T1, T2, T3>;
  arena_t<ret_type> ret
      = fma(value_of(arena_x), value_of(arena_y), value_of(arena_z));
  reverse_pass_callback(
      internal::fma_reverse_pass(arena_x, arena_y, arena_z, ret));
  return ret_type(ret);
}

}  // namespace math
}  // namespace stan
#endif
