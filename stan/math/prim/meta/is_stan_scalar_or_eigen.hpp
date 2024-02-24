#ifndef STAN_MATH_PRIM_META_IS_STAN_SCALAR_OR_EIGEN_HPP
#define STAN_MATH_PRIM_META_IS_STAN_SCALAR_OR_EIGEN_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/** \ingroup type_trait
 * Extends std::true_type if all the provided types are either a Stan Scalar
 * type or a type inheriting from `EigenBase`.
 */
template <typename T>
using is_stan_scalar_or_eigen
    = bool_constant<is_stan_scalar<std::decay_t<T>>::value
                    || is_eigen<std::decay_t<T>>::value>;

/*! \ingroup require_stan_scalar_real */
/*! \defgroup stan_scalar_or_eigen_types stan_scalar_or_eigen  */
/*! \addtogroup stan_scalar_or_eigen_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_stan_scalar_or_eigen */
/*! @tparam T the type to check */
template <typename T>
using require_stan_scalar_or_eigen_t
    = require_t<is_stan_scalar_or_eigen<std::decay_t<T>>>;
/*! @} */

}  // namespace stan
#endif
