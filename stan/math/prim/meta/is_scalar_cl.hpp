#ifndef STAN_MATH_PRIM_META_IS_SCALAR_CL_HPP
#define STAN_MATH_PRIM_META_IS_SCALAR_CL_HPP

#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {
namespace math {
namespace opencl {

template <typename T>
class ScalarCl;

}  // namespace opencl
}  // namespace math

namespace internal {
template <typename T>
struct is_scalar_cl_impl : std::false_type {};

template <typename T>
struct is_scalar_cl_impl<math::opencl::ScalarCl<T>> : std::true_type {
  using type = T;
};
}  // namespace internal

/** \ingroup type_traits
 * Checks if the decayed type of T is an `opencl::ScalarCl`, a scalar that
 * lives on the OpenCL device.
 */
template <typename T>
struct is_scalar_cl : internal::is_scalar_cl_impl<std::decay_t<T>> {};

namespace internal {
template <typename T>
struct is_prim_scalar_cl_impl : std::false_type {};

template <typename T>
struct is_prim_scalar_cl_impl<math::opencl::ScalarCl<T>>
    : std::is_arithmetic<T> {};

template <typename T>
struct is_rev_scalar_cl_impl : std::false_type {};

template <typename T>
struct is_rev_scalar_cl_impl<math::opencl::ScalarCl<T>> : is_var<T> {};
}  // namespace internal

/** \ingroup type_traits
 * Checks if the decayed type of T is an `opencl::ScalarCl` holding an
 * arithmetic value (`opencl::ScalarCl<double>`).
 */
template <typename T>
struct is_prim_scalar_cl : internal::is_prim_scalar_cl_impl<std::decay_t<T>> {};

/** \ingroup type_traits
 * Checks if the decayed type of T is an `opencl::ScalarCl` holding a var
 * (`opencl::ScalarCl<var>`).
 */
template <typename T>
struct is_rev_scalar_cl : internal::is_rev_scalar_cl_impl<std::decay_t<T>> {};

/** \ingroup type_traits
 * Specialization of `scalar_type` for `opencl::ScalarCl<T>`, which is `T`.
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_scalar_cl<T>::value>> {
  using type = typename internal::is_scalar_cl_impl<std::decay_t<T>>::type;
};

/** \ingroup type_traits
 * Specialization of `value_type` for `opencl::ScalarCl<T>`, which is `T`.
 */
template <typename T>
struct value_type<T, std::enable_if_t<is_scalar_cl<T>::value>> {
  using type = typename internal::is_scalar_cl_impl<std::decay_t<T>>::type;
};

/*! \ingroup matrix_cl_group */
/*! \defgroup scalar_cl_types scalar_cl  */
/*! \addtogroup scalar_cl_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_scalar_cl */
/*! @tparam T the type to check */
template <typename T>
using require_scalar_cl_t = require_t<is_scalar_cl<T>>;

/*! \brief Require type does not satisfy @ref is_scalar_cl */
/*! @tparam T the type to check */
template <typename T>
using require_not_scalar_cl_t = require_not_t<is_scalar_cl<T>>;

/*! \brief Require all of the types satisfy @ref is_scalar_cl */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_scalar_cl_t = require_all_t<is_scalar_cl<Types>...>;

/*! \brief Require any of the types satisfy @ref is_scalar_cl */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_scalar_cl_t = require_any_t<is_scalar_cl<Types>...>;

/*! \brief Require none of the types satisfy @ref is_scalar_cl */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_scalar_cl_t = require_all_not_t<is_scalar_cl<Types>...>;

/*! \brief Require type satisfies @ref is_prim_scalar_cl */
/*! @tparam T the type to check */
template <typename T>
using require_prim_scalar_cl_t = require_t<is_prim_scalar_cl<T>>;

/*! \brief Require type satisfies @ref is_rev_scalar_cl */
/*! @tparam T the type to check */
template <typename T>
using require_rev_scalar_cl_t = require_t<is_rev_scalar_cl<T>>;
/*! @} */

}  // namespace stan
#endif
