#ifndef STAN_MATH_REV_CORE_SAVE_VARIS_HPP
#define STAN_MATH_REV_CORE_SAVE_VARIS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>

#include <utility>
#include <vector>

namespace stan {
namespace math {

template <typename T, typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, const var_value<T>& x,
                                   Pargs&&... args);

template <typename T, typename VarVec,
          require_std_vector_vt<is_var, VarVec>* = nullptr, typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, VarVec&& x,
                                   Pargs&&... args);

template <typename T, typename VecContainer,
          require_std_vector_st<is_var, VecContainer>* = nullptr,
          require_std_vector_vt<is_container, VecContainer>* = nullptr,
          typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, VecContainer&& x,
                                   Pargs&&... args);

template <typename T, typename EigT, require_eigen_vt<is_var, EigT>* = nullptr,
          typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, EigT&& x,
                                   Pargs&&... args);

template <typename T, typename Arith,
          require_arithmetic_t<scalar_type_t<Arith>>* = nullptr,
          typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, Arith&& x,
                                   Pargs&&... args);
template <typename T>
inline vari_value<T>**& save_varis(vari_value<T>**& dest);

/**
 * End save_varis recursion and return pointer
 *
 * @param dest Pointer
 */
template <typename T>
inline vari_value<T>**& save_varis(vari_value<T>**& dest) {
  return dest;
}

/**
 * Save the vari pointer in x into the memory pointed to by dest,
 *   increment the dest storage pointer,
 *   recursively call save_varis on the rest of the arguments,
 *   and return the final value of the dest storage pointer.
 *
 * @tparam Pargs Types of remaining arguments
 * @param[in, out] dest Pointer to where vari pointers are saved
 * @param[in] x A var
 * @param[in] args Additional arguments to have their varis saved
 * @return Final position of dest pointer
 */
template <typename T, typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, const var_value<T>& x,
                                   Pargs&&... args) {
  *dest = x.vi_;
  return save_varis(dest + 1, std::forward<Pargs>(args)...);
}

/**
 * Save the vari pointers in x into the memory pointed to by dest,
 *   increment the dest storage pointer,
 *   recursively call save_varis on the rest of the arguments,
 *   and return the final value of the dest storage pointer.
 *
 * @tparam VarVec A variant of std::vector<var>
 * @tparam Pargs Types of remaining arguments
 * @param[in, out] dest Pointer to where vari pointers are saved
 * @param[in] x A std::vector of vars
 * @param[in] args Additional arguments to have their varis saved
 * @return Final position of dest pointer
 */
template <typename T, typename VarVec, require_std_vector_vt<is_var, VarVec>*,
          typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, VarVec&& x,
                                   Pargs&&... args) {
  for (int i = 0; i < x.size(); ++i) {
    dest[i] = x[i].vi_;
  }
  return save_varis(dest + x.size(), std::forward<Pargs>(args)...);
}

/**
 * Save the vari pointers in x into the memory pointed to by dest,
 *   increment the dest storage pointer,
 *   recursively call save_varis on the rest of the arguments,
 *   and return the final value of the dest storage pointer.
 *
 * @tparam VecContainer std::vector<T> where T is another type containing vars
 * @tparam Pargs Types of remaining arguments
 * @param[in, out] dest Pointer to where vari pointers are saved
 * @param[in] x A std::vector of containers containing of vars
 * @param[in] args Additional arguments to have their varis saved
 * @return Final position of dest pointer
 */
template <typename T, typename VecContainer,
          require_std_vector_st<is_var, VecContainer>*,
          require_std_vector_vt<is_container, VecContainer>*, typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, VecContainer&& x,
                                   Pargs&&... args) {
  for (size_t i = 0; i < x.size(); ++i) {
    dest = save_varis(dest, x[i]);
  }
  return save_varis(dest, std::forward<Pargs>(args)...);
}

/**
 * Save the vari pointers in x into the memory pointed to by dest,
 *   increment the dest storage pointer,
 *   recursively call save_varis on the rest of the arguments,
 *   and return the final value of the dest storage pointer.
 *
 * @tparam EigT An Eigen type with var value type
 * @tparam Pargs Types of remaining arguments
 * @param[in, out] dest Pointer to where vari pointers are saved
 * @param[in] x An Eigen container of vars
 * @param[in] args Additional arguments to have their varis saved
 * @return Final position of dest pointer
 */
template <typename T, typename EigT, require_eigen_vt<is_var, EigT>*,
          typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, EigT&& x,
                                   Pargs&&... args) {
  for (int i = 0; i < x.size(); ++i) {
    dest[i] = x(i).vi_;
  }
  return save_varis(dest + x.size(), std::forward<Pargs>(args)...);
}

/**
 * Ignore arithmetic types.
 *
 * Recursively call save_varis on the rest of the arguments
 *   and return the final value of the dest storage pointer.
 *
 * @tparam Arith An arithmetic type
 * @tparam Pargs Types of remaining arguments
 * @param[in, out] dest Pointer to where vari pointers are saved
 * @param[in] x An argument not containing vars
 * @param[in] args Additional arguments to have their varis saved
 * @return Final position of dest pointer
 */
template <typename T, typename Arith,
          require_arithmetic_t<scalar_type_t<Arith>>*, typename... Pargs>
inline vari_value<T>**& save_varis(vari_value<T>**& dest, Arith&& x,
                                   Pargs&&... args) {
  return save_varis(dest, std::forward<Pargs>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
