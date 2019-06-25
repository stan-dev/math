#ifndef STAN_MATH_PRIM_ARR_META_IS_CONSTANT_HPP
#define STAN_MATH_PRIM_ARR_META_IS_CONSTANT_HPP

#include <stan/math/prim/scal/meta/is_constant.hpp>
#include <vector>

namespace stan {
/**
 * Defines a public enum named value and sets it to true
 * if the type of the elements in the provided std::vector
 * is constant, false otherwise. This is used in
 * the is_constant_all metaprogram.
 * @tparam type of the elements in the std::vector
 */
template <typename T>
struct is_constant<std::vector<T> > {
  enum { value = is_constant<T>::value };
};

}  // namespace stan
#endif
