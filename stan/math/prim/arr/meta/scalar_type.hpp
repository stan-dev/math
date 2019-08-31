#ifndef STAN_MATH_PRIM_ARR_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_ARR_META_SCALAR_TYPE_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <vector>

namespace stan {
template <typename T>
struct scalar_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  typedef typename scalar_type<typename T::value_type>::type type;
};

}  // namespace stan
#endif
