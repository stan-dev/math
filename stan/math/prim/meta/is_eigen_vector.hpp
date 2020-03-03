#ifndef STAN_MATH_PRIM_META_IS_EIGEN_VECTOR_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_VECTOR_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_eigen_col_vector.hpp>
#include <stan/math/prim/meta/is_eigen_row_vector.hpp>
#include <type_traits>
#include <vector>


namespace stan {
  /** \ingroup type_trait
   * If the input type T is an eigen matrix with 1 column or 1 row at compile time
   * this has a static member with a value of true. Else this has a static
   * member with a value of false.
   */
  template <typename T>
  struct is_eigen_vector : bool_constant<is_eigen_col_vector<T>::value
                                         || is_eigen_row_vector<T>::value> {};

}

#endif
