#ifndef STAN_MATH_PRIM_SCAL_META_CONTAINS_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_CONTAINS_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/arr/meta/or.hpp>

namespace stan {

template <typename... T>
using contains_vector = or_<is_vector<T>...>;

}  // namespace stan
#endif
