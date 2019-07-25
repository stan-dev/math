#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_SAME_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_SAME_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T, typename S>
using enable_if_same = std::enable_if_t<std::is_same<T, S>::value>;

template <typename T, typename S>
using enable_if_not_same = std::enable_if_t<!std::is_same<T, S>::value>;

}  // namespace stan
#endif
