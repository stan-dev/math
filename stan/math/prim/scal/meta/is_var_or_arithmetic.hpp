#ifndef STAN_MATH_PRIM_SCAL_META_IS_VAR_OR_ARITHMETIC_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VAR_OR_ARITHMETIC_HPP

#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <type_traits>

namespace stan {

/**
 * Defines a public enum named value which is defined to be true (1)
 * if the type is either var or an aritmetic type
 * and false (0) otherwise.
 */
template <typename T>
struct is_var_or_arithmetic : std::integral_constant<bool,
 is_var<T>::value || std::is_arithmetic<T>::value>{};

// Helper Class to check if all input types are var or arithmetic
template <typename... T>
using is_all_var_or_arithmetic
    = math::conjunction<is_var_or_arithmetic<T>...>;

}  // namespace stan
#endif
