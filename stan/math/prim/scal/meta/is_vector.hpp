#ifndef STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP

#include <type_traits>

namespace stan {

template <typename T, typename = void>
struct is_vector : std::false_type {};

// We do treat pointers as vectors
template <typename T>
struct is_vector<T, std::enable_if_t<std::is_pointer<T>::value>> : std::true_type {};

template <typename T, typename = void>
struct is_std_vector : std::false_type {};

// Check whether type is an eigen col vector
template <typename T, typename = void>
struct is_eigen_col_vector : std::false_type {};

// Check whether type is an eigen row vector
template <typename T, typename = void>
struct is_eigen_row_vector : std::false_type {};

// Checks whether decayed type is an eigen vector
template <typename T, typename = void>
struct is_eigen_vector : std::false_type {};

}  // namespace stan
#endif
