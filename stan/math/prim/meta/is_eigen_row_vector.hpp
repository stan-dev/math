#ifndef STAN_MATH_PRIM_META_IS_EIGEN_ROW_VECTOR_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_ROW_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <type_traits>


namespace stan {
  namespace internal {
    /** \ingroup type_trait
     * Underlying implimentation for detecting if an Eigen Matrix is a row vector.
     */
    template <typename T, bool = stan::is_eigen<T>::value>
    struct is_eigen_row_vector_impl
        : std::integral_constant<bool, std::decay_t<T>::RowsAtCompileTime == 1> {};

    /** \ingroup type_trait
     * Specialization for when type is not an eigen vector.
     */
    template <typename T>
    struct is_eigen_row_vector_impl<T, false> : std::false_type {};

  }
  /** \ingroup type_trait
   * If the input type T is an eigen matrix with 1 column at compile time this
   * has a static member with a value of true. Else this has a static
   * member with a value of false.
   */
  template <typename T>
  struct is_eigen_row_vector : internal::is_eigen_row_vector_impl<T> {};

}

#endif
