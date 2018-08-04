#ifndef STAN_MATH_REV_CORE_BUILD_DOUBLE_ARRAY_HPP
#define STAN_MATH_REV_CORE_BUILD_DOUBLE_ARRAY_HPP

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
 * @return The newly populated array
 */
inline double* build_double_array(vari** array, int size) {
  double* double_array
      = ChainableStack::instance().memalloc_.alloc_array<double>(size);

  for (int i = 0; i < size; ++i)
    double_array[i] = array[i]->val_;

  return double_array;
}

/**
 * Make a copy of the input array allocated in the autodiff arena
 *
 * @param array Input array
 * @param size Number of doubles in input
 * @return Copy of input array
 */
inline double* build_double_array(const double* array, int size) {
  double* copy = ChainableStack::instance().memalloc_.alloc_array<double>(size);

  std::copy(array, array + size, copy);

  return copy;
}

}  // namespace math
}  // namespace stan
#endif
