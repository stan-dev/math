#ifndef STAN_MATH_PRIM_SCAL_META_SCALAR_TYPE_PRE_HPP
#define STAN_MATH_PRIM_SCAL_META_SCALAR_TYPE_PRE_HPP

#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>

namespace stan {
namespace internal {
template <bool is_vec, typename T, typename T_container>
struct scalar_type_helper_pre {
  typedef T_container type;
};

template <typename T, typename T_container>
struct scalar_type_helper_pre<true, T, T_container> {
  typedef
      typename scalar_type_helper_pre<is_vector<scalar_type_t<T>>::value,
                                      scalar_type_t<T>,
                                      scalar_type_t<T_container>>::type type;
};
}  // namespace internal

/**
 * Metaprogram structure to determine the type of first container of
 * the base scalar type of a template argument.
 *
 * @tparam T Type of object.
 */
template <typename T>
struct scalar_type_pre {
  typedef typename internal::scalar_type_helper_pre<
      is_vector<scalar_type_t<T>>::value, scalar_type_t<T>, T>::type type;
};

}  // namespace stan
#endif
