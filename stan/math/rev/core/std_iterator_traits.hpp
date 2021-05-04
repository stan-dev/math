#ifndef STAN_MATH_REV_CORE_STD_ITERATOR_TRAITS_HPP
#define STAN_MATH_REV_CORE_STD_ITERATOR_TRAITS_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>
#include <iterator>

namespace std {
/**
 * Specialization of iterator traits for Stan math.  These all take
 * the form of typedefs.
 */
template <typename Iter>
struct iterator_traits<stan::math::var_value<Iter>> {
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
  typedef stan::math::var_value<Iter> value_type;

  /**
   * Type of pointer to variables.
   */
  typedef stan::math::var_value<Iter>* pointer;

  /**
   * Type of reference to variables.
   */
  typedef stan::math::var_value<Iter>& reference;
};
}  // namespace std

#endif
