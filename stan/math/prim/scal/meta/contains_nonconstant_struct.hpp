#ifndef STAN_MATH_PRIM_SCAL_META_CONTAINS_NONCONSTANT_STRUCT_HPP
#define STAN_MATH_PRIM_SCAL_META_CONTAINS_NONCONSTANT_STRUCT_HPP

#include <stan/math/prim/scal/meta/is_constant.hpp>
#include <stan/math/prim/arr/meta/or.hpp>
#include <type_traits>
namespace stan {

template<typename... T>
  using contains_nonconstant_struct =  !or_<is_constant<T>...>;

}  // namespace stan
#endif
