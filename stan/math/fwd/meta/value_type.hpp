#ifndef STAN_MATH_FWD_META_VALUE_TYPE_HPP
#define STAN_MATH_FWD_META_VALUE_TYPE_HPP

#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/fwd/core/fvar.hpp>

#include <type_traits>

namespace stan {

template <typename T>
struct value_type<T, require_fvar_t<T>> {
  using type = typename std::decay_t<T>::Scalar;
};

}  // namespace stan
#endif
