#ifndef STAN_MATH_FWD_CORE_STD_ITERATOR_TRAITS_HPP
#define STAN_MATH_FWD_CORE_STD_ITERATOR_TRAITS_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/prim/meta.hpp>
#include <iterator>

namespace std {
/**
 * Specialization of iterator traits for Stan math.  These all take
 * the form of typedefs.
 */
template <typename T>
struct iterator_traits<stan::math::fvar<T>> {
  /**
   * Iterator category for traits.
   */
  typedef random_access_iterator_tag iterator_category;

  /**
   * Type for difference between pointers.
   */
  typedef ptrdiff_t difference_type;

  /**
   * Type for value of pointer to values.
   */
  typedef stan::math::fvar<T> value_type;

  /**
   * Type of pointer to variables.
   */
  typedef stan::math::fvar<T>* pointer;

  /**
   * Type of reference to variables.
   */
  typedef stan::math::fvar<T>& reference;
};
}  // namespace std

#endif
