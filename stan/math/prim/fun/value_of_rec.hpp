#ifndef STAN_MATH_PRIM_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_PRIM_FUN_VALUE_OF_REC_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/partially_forward_as_tuple.hpp>
#include <complex>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the value of the specified scalar argument
 * converted to a double value. For types with a \ref base_type
 *  of double this is itself. For types with a \ref base_type 
 *  of \ref var or \ref fvar this the `value_of_rec` of the 
 *  value type until a double is found.
 *
 * <p>See the <code>primitive_value</code> function to
 * extract values without casting to <code>double</code>.
 * @tparam A container or scalar type
 * @param x The type whose values will be extracted.
 * @return An object with a \ref base_type of double.
 */
template <typename T>
inline constexpr decltype(auto) value_of_rec(T&& x) {
  using val_t = std::decay_t<T>;
  // ints are cast to doubles, types with base double are passed along
  if constexpr (std::is_integral_v<val_t> || std::is_floating_point_v<val_t>) {
    return static_cast<double>(x);
  } else if constexpr (std::is_floating_point_v<base_type_t<val_t>>) {
    if constexpr (std::is_rvalue_reference_v<T&&>) {
      return plain_type_t<T>(std::forward<T>(x));
    } else {
      return x;
    }
  } else if constexpr (is_complex<val_t>::value) {
    return std::complex<double>{value_of_rec(x.real()), value_of_rec(x.imag())};
  } else if constexpr (is_std_vector_v<val_t>) {
    promote_scalar_t<double, val_t> ret;
    ret.reserve(x.size());
    for (auto&& x_i : x) {
      ret.push_back(value_of_rec(std::forward<decltype(x_i)>(x_i)));
    }
    return ret;
  } else if constexpr (is_eigen_v<val_t>) {
    return make_holder(
      [](auto& m) {
        return m.unaryExpr([](auto x_i) { return value_of_rec(x_i); });
      },
      std::forward<T>(x));
  } else if constexpr (is_tuple_v<val_t>) {
    return stan::math::apply(
      [](auto&&... args) {
        return partially_forward_as_tuple(
            value_of_rec(std::forward<decltype(args)>(args))...);
      },
      std::forward<T>(x));
  } else if constexpr (is_var_v<val_t>) {
    return x.vi_->val_;
  } else if constexpr (is_fvar<val_t>::value) {
    return value_of_rec(x.val());
  }
}

}  // namespace math
}  // namespace stan

#endif
