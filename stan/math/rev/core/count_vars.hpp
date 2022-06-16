#ifndef STAN_MATH_REV_CORE_COUNT_VARS_HPP
#define STAN_MATH_REV_CORE_COUNT_VARS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>

#include <utility>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename VecVar, require_std_vector_vt<is_var, VecVar>* = nullptr,
          typename... Pargs>
inline size_t count_vars_impl(size_t count, VecVar&& x, Pargs&&... args);

template <typename VecContainer,
          require_std_vector_st<is_var, VecContainer>* = nullptr,
          require_std_vector_vt<is_container, VecContainer>* = nullptr,
          typename... Pargs>
inline size_t count_vars_impl(size_t count, VecContainer&& x, Pargs&&... args);

template <typename EigT, require_eigen_vt<is_var, EigT>* = nullptr,
          typename... Pargs>
inline size_t count_vars_impl(size_t count, EigT&& x, Pargs&&... args);

template <typename... Pargs>
inline size_t count_vars_impl(size_t count, const var& x, Pargs&&... args);

template <typename Arith, require_arithmetic_t<scalar_type_t<Arith>>* = nullptr,
          typename... Pargs>
inline size_t count_vars_impl(size_t count, Arith& x, Pargs&&... args);

inline size_t count_vars_impl(size_t count);

/**
 * Count the number of vars in x (a std::vector of vars),
 *  add it to the running total,
 *  count the number of vars in the remaining arguments
 *  and return the result.
 *
 * @tparam VecVar type of standard container holding vars
 * @tparam Pargs Types of remaining arguments
 * @param[in] count The current count of the number of vars
 * @param[in] x A std::vector holding vars.
 * @param[in] args objects to be forwarded to recursive call of
 * `count_vars_impl`
 */
template <typename VecVar, require_std_vector_vt<is_var, VecVar>*,
          typename... Pargs>
inline size_t count_vars_impl(size_t count, VecVar&& x, Pargs&&... args) {
  return count_vars_impl(count + x.size(), std::forward<Pargs>(args)...);
}

/**
 * Count the number of vars in x (a std::vector holding other containers),
 *  add it to the running total,
 *  count the number of vars in the remaining arguments
 *  and return the result.
 *
 * @tparam VecContainer std::vector holding arguments which contain Vars
 * @tparam Pargs Types of remaining arguments
 * @param[in] count The current count of the number of vars
 * @param[in] x A vector holding containers of vars
 * @param[in] args objects to be forwarded to recursive call of
 * `count_vars_impl`
 */
template <typename VecContainer, require_std_vector_st<is_var, VecContainer>*,
          require_std_vector_vt<is_container, VecContainer>*, typename... Pargs>
inline size_t count_vars_impl(size_t count, VecContainer&& x, Pargs&&... args) {
  for (auto&& x_iter : x) {
    count = count_vars_impl(count, x_iter);
  }
  return count_vars_impl(count, std::forward<Pargs>(args)...);
}

/**
 * Count the number of vars in x (an eigen container),
 *  add it to the running total,
 *  count the number of vars in the remaining arguments
 *  and return the result.
 *
 * @tparam EigT A type derived from `EigenBase`
 * @tparam Pargs Types of remaining arguments
 * @param[in] count The current count of the number of vars
 * @param[in] x An Eigen container holding vars
 * @param[in] args objects to be forwarded to recursive call of
 * `count_vars_impl`
 */
template <typename EigT, require_eigen_vt<is_var, EigT>*, typename... Pargs>
inline size_t count_vars_impl(size_t count, EigT&& x, Pargs&&... args) {
  return count_vars_impl(count + x.size(), std::forward<Pargs>(args)...);
}

/**
 * Add one to the running total number of vars,
 *  count the number of vars in the remaining arguments
 *  and return the result.
 *
 * @tparam Pargs Types of remaining arguments
 * @param[in] count The current count of the number of vars
 * @param[in] x A var
 * @param[in] args objects to be forwarded to recursive call of
 * `count_vars_impl`
 */
template <typename... Pargs>
inline size_t count_vars_impl(size_t count, const var& x, Pargs&&... args) {
  return count_vars_impl(count + 1, std::forward<Pargs>(args)...);
}

/**
 * Arguments without vars contribute zero to the total number of vars.
 *
 * Return the running total number of vars plus the number of
 *   vars in the remaining arguments.
 *
 * @tparam Arith An object that is either arithmetic or holds arithmetic
 * types
 * @tparam Pargs Types of remaining arguments
 * @param[in] count The current count of the number of vars
 * @param[in] x An arithmetic value or container
 * @param[in] args objects to be forwarded to recursive call of
 * `count_vars_impl`
 */
template <typename Arith, require_arithmetic_t<scalar_type_t<Arith>>*,
          typename... Pargs>
inline size_t count_vars_impl(size_t count, Arith& x, Pargs&&... args) {
  return count_vars_impl(count, std::forward<Pargs>(args)...);
}

/**
 * End count_vars_impl recursion and return total number of counted vars
 */
inline size_t count_vars_impl(size_t count) { return count; }
}  // namespace internal

/**
 * Count the number of vars in the input argument list
 *
 * @tparam Pargs Types of input arguments
 * @return Number of vars in input
 */
template <typename... Pargs>
inline size_t count_vars(Pargs&&... args) {
  return internal::count_vars_impl(0, std::forward<Pargs>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
