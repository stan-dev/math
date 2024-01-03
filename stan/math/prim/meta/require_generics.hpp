#ifndef STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

// STAN_ADD_REQUIRE_BINARY(same, std::is_same, require_std);
template <typename T, typename S>
using require_same_t
    = require_t<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename S>
using require_not_same_t
    = require_not_t<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename... Types>
using require_all_same_t
    = require_all_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

template <typename T, typename... Types>
using require_any_same_t
    = require_any_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

template <typename T, typename... Types>
using require_all_not_same_t
    = require_all_not_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

template <typename T, typename... Types>
using require_any_not_same_t
    = require_any_not_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

template <typename T, typename... Types>
using require_any_not_convertible_t = require_any_not_t<
    std::is_convertible<std::decay_t<T>, std::decay_t<Types>>...>;

// STAN_ADD_REQUIRE_BINARY_INNER(same, std::is_same, require_std);
template <typename T, typename S>
using require_st_same = require_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                               scalar_type_t<std::decay_t<S>>>>;

template <typename T, typename S>
using require_not_st_same
    = require_not_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                 scalar_type_t<std::decay_t<S>>>>;

template <typename T, typename... Types>
using require_all_st_same
    = require_all_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                 scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_not_st_same
    = require_any_not_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                     scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename S>
using require_vt_same = require_t<
    std::is_same<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<S>>>>;

template <typename T, typename... Types>
using require_all_vt_same
    = require_all_t<std::is_same<value_type_t<std::decay_t<T>>,
                                 value_type_t<std::decay_t<Types>>>...>;

// STAN_ADD_REQUIRE_BINARY(convertible, std::is_convertible, require_std);
template <typename T, typename S>
using require_convertible_t
    = require_t<std::is_convertible<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename S>
using require_not_convertible_t
    = require_not_t<std::is_convertible<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename... Types>
using require_all_convertible_t = require_all_t<
    std::is_convertible<std::decay_t<T>, std::decay_t<Types>>...>;

// STAN_ADD_REQUIRE_BINARY(assignable, std::is_assignable, require_std);
template <typename T, typename S>
using require_assignable_t
    = require_t<std::is_assignable<std::decay_t<T>, std::decay_t<S>>>;

// STAN_ADD_REQUIRE_BINARY_INNER(assignable, std::is_assignable, require_std);
template <typename T, typename S>
using require_st_assignable
    = require_t<std::is_assignable<scalar_type_t<std::decay_t<T>>,
                                   scalar_type_t<std::decay_t<S>>>>;

template <typename T, typename S>
using require_not_st_assignable
    = require_not_t<std::is_assignable<scalar_type_t<std::decay_t<T>>,
                                       scalar_type_t<std::decay_t<S>>>>;

template <typename T, typename... Types>
using require_all_st_assignable
    = require_all_t<std::is_assignable<scalar_type_t<std::decay_t<T>>,
                                       scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_st_assignable
    = require_any_t<std::is_assignable<scalar_type_t<std::decay_t<T>>,
                                       scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_all_not_st_assignable = require_all_not_t<std::is_assignable<
    scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_not_st_assignable = require_any_not_t<std::is_assignable<
    scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename S>
using require_vt_assignable
    = require_t<std::is_assignable<value_type_t<std::decay_t<T>>,
                                   value_type_t<std::decay_t<S>>>>;

template <typename T, typename S>
using require_not_vt_assignable
    = require_not_t<std::is_assignable<value_type_t<std::decay_t<T>>,
                                       value_type_t<std::decay_t<S>>>>;

template <typename T, typename... Types>
using require_all_vt_assignable
    = require_all_t<std::is_assignable<value_type_t<std::decay_t<T>>,
                                       value_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_vt_assignable
    = require_any_t<std::is_assignable<value_type_t<std::decay_t<T>>,
                                       value_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_all_not_vt_assignable = require_all_not_t<std::is_assignable<
    value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_not_vt_assignable = require_any_not_t<std::is_assignable<
    value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>;

// STAN_ADD_REQUIRE_BINARY(constructible, std::is_constructible, require_std);
template <typename T, typename S>
using require_constructible_t
    = require_t<std::is_constructible<std::decay_t<T>, std::decay_t<S>>>;

// STAN_ADD_REQUIRE_BINARY_INNER(constructible, std::is_constructible,
// require_std);
template <typename T, typename S>
using require_st_constructible
    = require_t<std::is_constructible<scalar_type_t<std::decay_t<T>>,
                                      scalar_type_t<std::decay_t<S>>>>;

template <typename T, typename S>
using require_not_st_constructible
    = require_not_t<std::is_constructible<scalar_type_t<std::decay_t<T>>,
                                          scalar_type_t<std::decay_t<S>>>>;

template <typename T, typename... Types>
using require_all_st_constructible = require_all_t<std::is_constructible<
    scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_st_constructible = require_any_t<std::is_constructible<
    scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_all_not_st_constructible
    = require_all_not_t<std::is_constructible<
        scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_not_st_constructible
    = require_any_not_t<std::is_constructible<
        scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename S>
using require_vt_constructible
    = require_t<std::is_constructible<value_type_t<std::decay_t<T>>,
                                      value_type_t<std::decay_t<S>>>>;

template <typename T, typename S>
using require_not_vt_constructible
    = require_not_t<std::is_constructible<value_type_t<std::decay_t<T>>,
                                          value_type_t<std::decay_t<S>>>>;

template <typename T, typename... Types>
using require_all_vt_constructible = require_all_t<std::is_constructible<
    value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_vt_constructible = require_any_t<std::is_constructible<
    value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_all_not_vt_constructible
    = require_all_not_t<std::is_constructible<
        value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>;

template <typename T, typename... Types>
using require_any_not_vt_constructible
    = require_any_not_t<std::is_constructible<
        value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>;

/*! \ingroup require_stan_scalar_real */
/*! \defgroup arithmetic_types arithmetic  */
/*! \addtogroup arithmetic_types */
/*! @{ */

/*! \brief Require type satisfies std::is_arithmetic */
template <typename T>
using require_arithmetic_t = require_t<std::is_arithmetic<std::decay_t<T>>>;

/*! \brief Require type does not satisfy std::is_arithmetic */
template <typename T>
using require_not_arithmetic_t
    = require_not_t<std::is_arithmetic<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy std::is_arithmetic */
template <typename... Types>
using require_all_arithmetic_t
    = require_all_t<std::is_arithmetic<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy std::is_arithmetic */
template <typename... Types>
using require_any_arithmetic_t
    = require_any_t<std::is_arithmetic<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy std::is_arithmetic */
template <typename... Types>
using require_all_not_arithmetic_t
    = require_all_not_t<std::is_arithmetic<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy std::is_arithmetic */
template <typename... Types>
using require_any_not_arithmetic_t
    = require_any_not_t<std::is_arithmetic<std::decay_t<Types>>...>;
/*! @} */


/*! \ingroup require_stan_scalar_real */
/*! \addtogroup arithmetic_types */
/*! @{ */

/*! \brief Require all of the value types satisfy std::is_arithmetic */
template <typename... Types>
using require_all_vt_arithmetic
    = require_all_t<std::is_arithmetic<value_type_t<std::decay_t<Types>>>...>;

/*! \brief Require at least one of the value types do not satisfy std::is_arithmetic */
template <typename... Types>
using require_any_not_vt_arithmetic = require_any_not_t<
    std::is_arithmetic<value_type_t<std::decay_t<Types>>>...>;

/*! \brief Require scalar type satisfies std::is_arithmetic */
template <typename T>
using require_st_arithmetic
    = require_t<std::is_arithmetic<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require scalar type does not satisfy std::is_arithmetic */
template <typename T>
using require_not_st_arithmetic
    = require_not_t<std::is_arithmetic<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require all of the scalar types satisfy std::is_arithmetic */
template <typename... Types>
using require_all_st_arithmetic
    = require_all_t<std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>;

/*! \brief Require any of the scalar types satisfy std::is_arithmetic */
template <typename... Types>
using require_any_st_arithmetic
    = require_any_t<std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>;

/*! \brief Any of the scalar types do not satisfy std::is_arithmetic */
template <typename... Types>
using require_any_not_st_arithmetic = require_any_not_t<
    std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>;
/*! @} */

/*! \ingroup require_stan_scalar_real */
/*! \defgroup floating_point_types floating_point  */
/*! \addtogroup floating_point_types */
/*! @{ */

/*! \brief Require type satisfies std::is_floating_point */
template <typename T>
using require_floating_point_t
    = require_t<std::is_floating_point<std::decay_t<T>>>;
/*! @} */

/*! \ingroup require_stan_scalar_real */
/*! \defgroup integral_types integral  */
/*! \addtogroup integral_types */
/*! @{ */

/*! \brief Require type satisfies std::is_integral */
template <typename T>
using require_integral_t = require_t<std::is_integral<std::decay_t<T>>>;
/*! @} */

/*! \ingroup require_stan_scalar_real */
/*! \addtogroup integral_types */
/*! @{ */

/*! \brief Require value type satisfies std::is_integral */
template <typename T>
using require_vt_integral
    = require_t<std::is_integral<value_type_t<std::decay_t<T>>>>;

/*! \brief Require scalar type satisfies std::is_integral */
template <typename T>
using require_st_integral
    = require_t<std::is_integral<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require scalar type does not satisfy std::is_integral */
template <typename T>
using require_not_st_integral
    = require_not_t<std::is_integral<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require none of the scalar types satisfy std::is_integral */
template <typename... Types>
using require_all_not_st_integral = require_all_not_t<
    std::is_integral<scalar_type_t<std::decay_t<Types>>>...>;
/*! @} */

}  // namespace stan
#endif
