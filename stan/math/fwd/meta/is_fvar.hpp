#ifndef STAN_MATH_FWD_META_IS_FVAR_HPP
#define STAN_MATH_FWD_META_IS_FVAR_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <type_traits>

namespace stan {

namespace internal {
template <typename T>
struct is_fvar_impl : std::false_type {};

template <typename T>
struct is_fvar_impl<math::fvar<T>> : std::true_type {};

}  // namespace internal

template <typename T>
struct is_fvar<T,
               std::enable_if_t<internal::is_fvar_impl<std::decay_t<T>>::value>>
    : std::true_type {};

}  // namespace stan
#endif
