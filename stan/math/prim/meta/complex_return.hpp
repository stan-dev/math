#ifndef STAN_MATH_PRIM_META_COMPLEX_RETURN_HPP
#define STAN_MATH_PRIM_META_COMPLEX_RETURN_HPP

#include <stan/math/prim/meta/promote_args.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>

namespace stan {
namespace math {

namespace internal {

template <typename U, typename V>
struct complex_return_binary {};

template <typename U, typename V>
struct complex_return_binary<std::complex<U>, std::complex<V>> {
  using type = std::complex<return_type_t<U, V>>;
};

template <typename U, typename V>
struct complex_return_binary<std::complex<U>, V> {
  using type = std::complex<return_type_t<U, V>>;
};

template <typename U, typename V>
struct complex_return_binary<U, std::complex<V>> {
  using type = std::complex<return_type_t<U, V>>;
};

template <typename... Args>
using complex_return_binary_t = typename complex_return_binary<Args...>::type;

}  // namespace internal

/** \ingroup type_trait
 * Template metaprogram determining type of complex return for
 * a function of one or more arguments.  All returns will be of
 * minimal type `std::complex<double>` even for all `double`
 * arguments.  Otherwise, template arguments to `std::complex` will be
 * combined using the metaprogram `stan::math::return_type_t`.
 *
 * @tparam T type of first argument
 * @tparam Ts types of remaining arguments
 */
template <typename T, typename... Ts>
struct complex_return {
  using type
      = internal::complex_return_binary_t<T,
                                          typename complex_return<Ts...>::type>;
};

template <typename T>
struct complex_return<T> {
  using type = internal::complex_return_binary_t<T, std::complex<double>>;
};

template <typename... Args>
using complex_return_t = typename complex_return<Args...>::type;

}  // namespace math
}  // namespace stan

#endif
