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
    T, require_all_t<is_var<scalar_type_t<T>>,
                     math::disjunction<is_eigen<T>, is_eigen<value_type_t<T>>>>>
    : std::true_type {};

}  // namespace stan
#endif
