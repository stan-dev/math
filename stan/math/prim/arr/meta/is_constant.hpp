#ifndef STAN_MATH_PRIM_ARR_META_IS_CONSTANT_HPP
#define STAN_MATH_PRIM_ARR_META_IS_CONSTANT_HPP

#include <stan/math/prim/scal/meta/is_constant.hpp>
#include <vector>

namespace stan {
/**
 * Defines a public enum named value and sets it to true(1)
 * if the type of the elements in the provided std::vector
 * is a constant struct, false(0) otherwise. This helper
 * struct is used in the is_constant_struct metaprogram.
 * @tparam type of the elements in the std::vector
 */
template <typename T>
struct is_constant_all<std::vector<T> > {
  enum { value = is_constant_all<T>::value };
};

}  // namespace stan
#endif
