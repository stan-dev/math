#ifndef STAN_MATH_REV_CORE_BUILD_DOUBLE_ARRAY_IF_NECESSARY_HPP
#define STAN_MATH_REV_CORE_BUILD_DOUBLE_ARRAY_IF_NECESSARY_HPP

#include <stan/math/rev/core/var.hpp>
#include <limits>

namespace stan {
namespace math {
/**
 * Allocate space for and populate an array in memory with the values
 * of the vars pointed to by array. Return the allocated array
 *
 * @param array Array of vars to extract values from
 * @param size Number of vars in input
 * @return The new allocated array
 */
double* build_double_array_if_necessary(var* array, int size) {
  double* double_array
      = ChainableStack::instance().memalloc_.alloc_array<double>(size);

  for (int i = 0; i < size; ++i)
    double_array[i] = array[i].val();

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
double* build_double_array_if_necessary(double* array, int size) {
  return array;
}
}  // namespace math
}  // namespace stan
#endif
