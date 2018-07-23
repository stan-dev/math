#ifndef STAN_MATH_REV_CORE_BUILD_VARI_POINTER_ARRAY_IF_NECESSARY_HPP
#define STAN_MATH_REV_CORE_BUILD_VARI_POINTER_ARRAY_IF_NECESSARY_HPP

#include <stan/math/rev/core/var.hpp>
#include <limits>

namespace stan {
namespace math {
/**
 * Allocate space for and populate an array in memory with the values
 * of the varis pointed to by the elements of array. Return the allocated array
 *
 * @param array Array of pointers to varis to extract values from
 * @param size Number of elements in input
 * @return The new allocated array
 */
inline vari** build_vari_pointer_array_if_necessary(const var* array,
                                                    int size) {
  vari** vari_pointer_array
      = ChainableStack::instance().memalloc_.alloc_array<vari*>(size);

  for (int i = 0; i < size; ++i)
    vari_pointer_array[i] = array[i].vi_;

  return vari_pointer_array;
}

/**
 * build_vari_pointer_array_if_necessary just returns the input if it is an
 * array of doubles
 *
 * @param array Input array
 * @param size Number of doubles in input
 * @return The input
 */
inline const double* build_vari_pointer_array_if_necessary(const double* array,
                                                           int size) {
  return array;
}
}  // namespace math
}  // namespace stan
#endif
