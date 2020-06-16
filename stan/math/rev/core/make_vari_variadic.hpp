#ifndef STAN_MATH_REV_CORE_MAKE_VARI_VARIANT_HPP
#define STAN_MATH_REV_CORE_MAKE_VARI_VARIANT_HPP

#include <stan/math/rev/meta.hpp>
#include <boost/variant2/variant.hpp>
#include <tuple>

namespace stan {
namespace math {
template <typename T, typename = void>
class vari_value;

namespace internal {
// Make one tuple out of N tuples
template <typename... Types>
struct concat_tuple {
  using type = decltype(std::tuple_cat(std::declval<Types>()...));
};
// helper
template <typename... Types>
using concat_tuple_t = typename concat_tuple<Types...>::type;

// All the vari_value Matrix types
template <typename T>
struct vari_eigen_mat_tuple {
    using type = std::tuple<
     vari_value<Eigen::Matrix<T, -1, -1>>,
     vari_value<Eigen::Matrix<T, 1, -1>>,
     vari_value<Eigen::Matrix<T, -1, 1>>,
     vari_value<Eigen::Matrix<T, 1, 1>>>;
};

template <typename T>
using vari_eigen_mat_tuple_t = typename vari_eigen_mat_tuple<T>::type;

template <typename T>
struct vari_eigen_arr_tuple {
    using type = std::tuple<
     vari_value<Eigen::Array<T, -1, -1>>,
     vari_value<Eigen::Array<T, 1, -1>>,
     vari_value<Eigen::Array<T, -1, 1>>,
     vari_value<Eigen::Array<T, 1, 1>>>;
};

template <typename T>
using vari_eigen_arr_tuple_t = typename vari_eigen_arr_tuple<T>::type;

template <typename T>
struct vari_eigen_sparsemat_tuple {
    using type = std::tuple<
     vari_value<Eigen::SparseMatrix<T>>,
     vari_value<Eigen::SparseVector<T>>>;
};

template <typename T>
using vari_eigen_sparsemat_tuple_t = typename vari_eigen_sparsemat_tuple<T>::type;

// All of the vari_value mixes
template <typename T>
struct vari_combinations {
    using type = concat_tuple_t<
      vari_eigen_mat_tuple_t<T>,
      vari_eigen_arr_tuple_t<T>,
//      vari_eigen_sparsemat_tuple_t<T>,
      std::tuple<vari_value<T>>>;
};

template <typename T>
using vari_combinations_t = typename vari_combinations<T>::type;


template <typename T = void, typename... Types>
struct make_vari_types;

 // break condition on default type which ends the recursion
template <>
struct make_vari_types<void> {
     using type = std::tuple<>;
};

// alias template
template <typename... Types>
using make_vari_types_t = typename make_vari_types<Types...>::type;

// Take all the base types T and make their vari_value equivalents
template <typename T, typename... Types>
struct make_vari_types {
  using type = concat_tuple_t<vari_combinations_t<T>, make_vari_types_t<Types...>>;
};

template <typename... Types>
struct make_vari_variadic {};

// Take the tuple of all the types and place them in boost::variant
template <typename... Types>
struct make_vari_variadic<std::tuple<Types...>> {
    using type = boost::variant2::variant<Types*...>;
};
}
// Create boost variant holding pointers to all the vari_value types
template <typename... Types>
using make_vari_variadic_t = typename internal::make_vari_variadic<internal::make_vari_types_t<Types...>>::type;

using vari_variant = make_vari_variadic_t<float, long double, double>;
}
}
#endif
