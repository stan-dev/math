#ifndef STAN_MATH_REV_META_CONDITIONAL_SEQUENCE_HPP
#define STAN_MATH_REV_META_CONDITIONAL_SEQUENCE_HPP

#include <utility>
#include <stddef.h>

namespace stan {

/**
 * End the recursion for making a compile time integer sequence
 * @tparam J a size_t valued nttp that is tossed.
 * @tparam I A parameter pack of size_t nttps that fill out the returned
 *  index sequence.
 */
template <template <typename...> class Checker, typename T, size_t J,
          size_t... I>
constexpr auto conditional_sequence(Checker<T> x,
                                    std::index_sequence<J, I...> /* ignore */) {
  return std::index_sequence<I...>{};
}

/**
 * Create a conditinal compile time integer sequence that increments
 * based on the number of Eigen types with numeric scalars in the
 * arguments.
 * @tparam J The starting value for the integer sequence.
 * @tparam I A parameter pack representing the integer sequence that is passed
 *  recursivly through calls to this function.
 * @tparam T The first input type peeled off through recursive calls. If this
 * type is an Eigen matrix with an arithmetic scalar the next recursive call
 * gets a +1 in the integer sequence.
 * @tparam Types The remaining parameter pack of types used to calculate the
 *   the integer sequence.
 *
 * For a given call such as
 * ```
 * make_cond_sequence(std::integer_sequence<0>, var, var,
 *   Matrix<double>, var, Matrix<double>, Matrix<double>, var)
 * ```
 * This returns the output
 * `index_sequence<0, 0, 0, 1, 1, 2, 2>`
 * which can be used by `make_op_vari` to extract the correct pointer to memory
 * in the stack for each eigen matrix with arithmetic scalars.
 */
template <template <typename...> class Checker, typename T1, size_t J,
          size_t... I, typename T, typename... Types>
constexpr auto conditional_sequence(Checker<T1> toss,
                                    std::index_sequence<J, I...> /* ignore */,
                                    T&& x, Types&&... args) {
  return conditional_sequence(
      toss, std::index_sequence<J + Checker<std::decay_t<T>>::value, I..., J>{},
      args...);
}
}  // namespace stan
#endif
