#ifndef STAN_MATH_REV_CORE_BUILD_DOUBLE_ARRAY_IF_NECESSARY_HPP
#define STAN_MATH_REV_CORE_BUILD_DOUBLE_ARRAY_IF_NECESSARY_HPP

#include <stan/math/rev/core/vari.hpp>
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
inline const double* build_double_array_if_necessary(vari** array, int size) {
  double* double_array
      = ChainableStack::instance().memalloc_.alloc_array<double>(size);

  for (int i = 0; i < size; ++i)
    double_array[i] = array[i]->val_;

  return double_array;
}

/**
 * build_double_array_if_necessary does nothing if the input is already an array
 * of doubles. Return the input
 *
 * @param array Input array
 * @param size Number of doubles in input
 * @return The input
 */
inline const double* build_double_array_if_necessary(const double* array,
                                                     int size) {
  return array;
}
}  // namespace math
}  // namespace stan
#endif
