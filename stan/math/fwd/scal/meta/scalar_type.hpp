#ifndef STAN_MATH_FWD_SCAL_META_SCALAR_TYPE_HPP
#define STAN_MATH_FWD_SCAL_META_SCALAR_TYPE_HPP

#include <stan/math/fwd/scal/meta/is_fvar.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>
template <typename T>
class fvar;
namespace stan {

/**
 * Specialization of scalar_type for fvar
 */
template <template <class> class K, typename T>
struct scalar_type<K<T>, std::enable_if_t<is_fvar<K<T>>::value>> {
  typedef fvar<T> type;
};

}  // namespace stan
#endif
