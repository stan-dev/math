#ifndef STAN_TEST_UNIT_MATH_MIX_UTIL_HPP
#define STAN_TEST_UNIT_MATH_MIX_UTIL_HPP

#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/rev/fun/util.hpp>

namespace stan::math::test {
/**
 * For autodiff types, return the double value of the input's value
 * @tparam T Type with a `auto val()` member function or a nonad scalar type
 * @param x Scalar value
 * @return The double value of the input's value
 */
template <typename T>
inline double get_val(T&& x) {
  if constexpr (stan::is_fvar_v<T> || stan::is_var_v<T>) {
    return get_val(x.val());
  } else if constexpr (stan::is_stan_scalar_v<T>) {
    return x;
  } else {
    static_assert(sizeof(T*) == 0, "This function only accepts scalar types!");
  }
}
/**
 * Recurse through tuples and containers and apply a function to scalars.
 * @tparam F A functor with operator() that takes scalar types and returns void.
 * @tparam Types A variadic pack of types.
 * @param f The functor to apply to scalars.
 * @param x args to apply f to.
 */
template <typename F, typename... Types>
inline void recursive_for_each(F&& f, Types&&... args) {
  if constexpr (std::conjunction_v<stan::math::is_tuple<Types>...>) {
    stan::math::for_each(
        [&f](auto&&... args_i) { recursive_for_each(f, args_i...); }, args...);
  } else {
    if constexpr (std::conjunction_v<stan::is_std_vector<Types>...>) {
      const auto max_size = stan::math::max_size(args...);
      for (Eigen::Index i = 0; i < max_size; ++i) {
        if constexpr (std::conjunction_v<
                          stan::is_stan_scalar<value_type_t<Types>>...>) {
          f(args[i]...);
        } else {
          recursive_for_each(f, args[i]...);
        }
      }
    } else if constexpr (std::conjunction_v<stan::is_eigen<Types>...>) {
      const auto max_size = stan::math::max_size(args...);
      for (Eigen::Index i = 0; i < max_size; ++i) {
        f(args(i)...);
      }
    } else if constexpr (std::conjunction_v<stan::is_stan_scalar<Types>...>) {
      f(args...);
    }
  }
}

}  // namespace stan::math::test

#endif
