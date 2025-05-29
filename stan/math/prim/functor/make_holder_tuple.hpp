#ifndef STAN_MATH_PRIM_FUNCTOR_MAKE_HOLDER_TUPLE_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAKE_HOLDER_TUPLE_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * @brief Helper template to deduce the correct type for tuple elements.
 *
 * This template determines the type to store in a tuple when forwarding
 * arguments.
 * - If the argument type T is an rvalue reference, the resulting type is the
 * decayed type (i.e. cv-qualified types and references are removed) so that the
 * tuple element becomes a value.
 * - If the argument type T is not an rvalue reference (typically deduced as an
 * lvalue reference), then T&& will collapse to an lvalue reference, preserving
 * the reference.
 *
 * @tparam T The type to deduce.
 */
template <typename T>
struct deduce_cvr {
  using type
      = std::conditional_t<std::is_rvalue_reference_v<T>, std::decay_t<T>, T&&>;
};

template <typename T>
using deduce_cvr_t = typename deduce_cvr<T>::type;
}  // namespace internal
/**
 * Holds ownership of rvalues and forwards lvalues into a tuple.
 *
 * Constructs a tuple from the provided arguments such that:
 * - If an argument is an xvalue (an rvalue), the tuple element will be a
 * decayed value.
 * - If an argument is an lvalue (or const lvalue), the tuple element will
 * maintain its reference type.
 *
 * This behavior ensures that temporaries are stored by value in the tuple while
 * lvalues are preserved as references. It is similar in intent to the `Holder`
 * class in behavior, but for tuples instead of Eigen types.
 * It is the opposite of `std::forward_as_tuple`, with the difference in
 * handling rvalues. `std::forward_as_tuple` does not extend object lifetimes,
 * so when an rvalue is passed to `std::forward_as_tuple`, the resulting tuple
 * element will be a reference to a temporary that is destroyed at the end of
 * the statement. This function ensures that rvalues are stored by value in the
 * tuple, extending their lifetimes.
 *
 * @tparam Types Parameter pack representing the types of the arguments.
 * @param args The arguments to forward into the tuple.
 * @return A tuple containing the forwarded arguments with types deduced via
 * deduce_cvr_t.
 *
 * @note The function is declared constexpr. It is noexcept when, for each input
 * type, the type value coming in is either an lvalue reference or an rvalue
 * type whose move constructor is nothrow
 */
template <typename... Types>
inline constexpr auto make_holder_tuple(Types&&... args) {
  if constexpr (sizeof...(Types) == 0) {
    return std::tuple<>{};
  } else {
    return std::tuple<internal::deduce_cvr_t<Types&&>...>{
        std::forward<Types>(args)...};
  }
}
}  // namespace math
}  // namespace stan
#endif
