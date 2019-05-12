#ifndef STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP

#include <stan/math/prim/scal/meta/disjunction.hpp>

namespace stan {

// FIXME: use boost::type_traits::remove_all_extents to
//        extend to array/ptr types

template <typename T>
struct is_vector_helper {
  enum { value = 0 };
  typedef T type;
};

template <typename... T>
using is_vector = math::disjunction<is_vector_helper<T>...>;
}  // namespace stan
#endif
