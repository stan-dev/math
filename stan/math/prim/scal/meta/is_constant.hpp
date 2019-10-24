#ifndef STAN_MATH_PRIM_SCAL_META_IS_CONSTANT_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_CONSTANT_HPP

#include <stan/math/prim/scal/meta/bool_constant.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <type_traits>

namespace stan {

/**
 * Metaprogramming struct to detect whether a given type is constant
 * in the mathematical sense (not the C++ <code>const</code>
 * sense). If the parameter type is constant, <code>value</code>
 * will be equal to <code>true</code>.
 *
 * The baseline implementation in this abstract base class is to
 * classify a type <code>T</code> as constant if it can be converted
 * (i.e., assigned) to a <code>double</code>.  This baseline should
 * be overridden for any type that should be treated as a variable.
 *
 * @tparam T Type being tested.
 */
template <typename T, typename = void>
struct is_constant : bool_constant<std::is_convertible<T, double>::value> {};

/**
 * Metaprogram defining an enum <code>value</code> which
 * is <code>true</code> if all of the type parameters
 * are constant (i.e., primtive types) and
 * <code>false</code> otherwise.
 */
template <typename... T>
using is_constant_all = math::conjunction<is_constant<T>...>;

}  // namespace stan
#endif
