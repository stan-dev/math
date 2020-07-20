#ifndef STAN_MATH_FWD_META_PARTIALS_TYPE_HPP
#define STAN_MATH_FWD_META_PARTIALS_TYPE_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta/is_fvar.hpp>
#include <stan/math/prim/meta/partials_type.hpp>
#include <type_traits>

namespace stan {

template <typename T>
struct partials_type<T, require_fvar_t<T>> {
  using type = typename std::decay_t<T>::Scalar;
};

}  // namespace stan
#endif
