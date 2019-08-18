#ifndef STAN_MATH_FWD_SCAL_META_SCALAR_TYPE_HPP
#define STAN_MATH_FWD_SCAL_META_SCALAR_TYPE_HPP

#include <stan/math/fwd/scal/meta/is_fvar.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>
namespace stan {

/**
 * Specialization of scalar_type for fvar
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_fvar<T>::value>> {
  typedef T type;
};

}  // namespace stan
#endif
