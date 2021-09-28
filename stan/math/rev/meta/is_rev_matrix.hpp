#ifndef STAN_MATH_REV_META_IS_REV_MATRIX_HPP
#define STAN_MATH_REV_META_IS_REV_MATRIX_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_rev_matrix.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either a type derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`
 */
template <typename T>
struct is_rev_matrix<
    T,
    require_all_t<is_var<scalar_type_t<T>>,
                  math::disjunction<
                      math::conjunction<is_var<T>, is_eigen<value_type_t<T>>>,
                      is_eigen<T>>>> : std::true_type {};

/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either a type derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of columns equal to 1.
 */
template <typename T>
struct is_rev_col_vector<
    T, require_all_t<is_var<scalar_type_t<T>>,
                     math::disjunction<is_eigen_col_vector<T>,
                                       is_eigen_col_vector<value_type_t<T>>>>>
    : std::true_type {};

/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either a type derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of rows equal to 1.
 */
template <typename T>
struct is_rev_row_vector<
    T, require_all_t<is_var<scalar_type_t<T>>,
                     math::disjunction<is_eigen_row_vector<T>,
                                       is_eigen_row_vector<value_type_t<T>>>>>
    : std::true_type {};

/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is either a type derived from `Eigen::EigenBase` with a `Scalar`
 *  type of `var_value<double>` or a `var_value<T>` where T is derived from
 * `Eigen::EigenBase`. And the type must have a compile time constant number
 *  of rows equal to 1.
 */
template <typename T>
struct is_rev_vector<T,
                     require_any_t<is_rev_col_vector<T>, is_rev_row_vector<T>>>
    : std::true_type {};

}  // namespace stan
#endif
