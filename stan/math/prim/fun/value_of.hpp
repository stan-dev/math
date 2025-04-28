#ifndef STAN_MATH_PRIM_FUN_VALUE_OF_HPP
#define STAN_MATH_PRIM_FUN_VALUE_OF_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/partially_forward_as_tuple.hpp>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the value of the specified argument.
 *  For types with a \ref base_type of double or int returns itself.
 *  For types with a \ref base_type of \ref var or \ref fvar
 *  the `value` member of their class is returned.
 *
 *  So for `std::complex<fvar<var>>` this will return
 *  a `std::complex<var>`. And for `std::vector<var>`
 *  this will return a `std:vector<double>` whose
 *  values are the `val_` members of the `var`s.
 *
 * <p>See the <code>primitive_value</code> function to
 * extract values without casting to <code>double</code>.
 * @tparam T A container or scalar type
 * @param x The object whose values will be extracted.
 * @return An object whose \ref scalar_type
 */
template <typename T>
inline constexpr decltype(auto) value_of(T&& x) {
  using val_t = std::decay_t<T>;
  // ints are cast to doubles, types with base double are passed along
  if constexpr (is_tuple_v<val_t>) {
    return stan::math::apply(
        [](auto&&... args) {
          return partially_forward_as_tuple(
              value_of(std::forward<decltype(args)>(args))...);
        },
        std::forward<T>(x));
  } else {
    constexpr bool is_base_float_or_int = std::is_floating_point_v<base_type_t<val_t>> || std::is_integral_v<base_type_t<val_t>>;
    if constexpr (std::is_integral_v<val_t> || std::is_floating_point_v<val_t>) {
      return val_t{x};
    } else if constexpr (is_base_float_or_int) {
      if constexpr (std::is_rvalue_reference_v<T&&>) {
        return plain_type_t<T>(std::forward<T>(x));
      } else {
        return x;
      }
    } else if constexpr (is_complex<val_t>::value) {
      return std::complex<double>{value_of(x.real()), value_of(x.imag())};
    } else if constexpr (is_std_vector_v<val_t>) {
      std::vector<
          plain_type_t<decltype(value_of(std::declval<value_type_t<T>>()))>>
          ret;
      ret.reserve(x.size());
      for (auto&& x_i : x) {
        ret.push_back(value_of(std::forward<decltype(x_i)>(x_i)));
      }
      return ret;
    } else if constexpr (is_eigen_v<val_t>) {
      return make_holder(
          [](auto& m) {
            return m.unaryExpr([](auto x_i) { return value_of(x_i); });
          },
          std::forward<T>(x));
    } else if constexpr (is_var_v<val_t>) {
      return x.vi_->val_;
    } else if constexpr (is_fvar<val_t>::value) {
      return x.val();
    }
  
  }
}

}  // namespace math
}  // namespace stan

#endif
