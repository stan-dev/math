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
inline vari** save_varis(vari** dest, const var& x) {
  *dest = x.vi_;
  return dest + 1;
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
template <typename VarVec, require_std_vector_vt<is_var, VarVec>* = nullptr>
inline vari** save_varis(vari** dest, VarVec&& x) {
  for (int i = 0; i < x.size(); ++i) {
    dest[i] = x[i].vi_;
  }
  return dest + x.size();
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
template <typename VecContainer, require_std_vector_st<is_var, VecContainer>* = nullptr,
          require_std_vector_vt<is_container, VecContainer>* = nullptr>
inline vari** save_varis(vari** dest, VecContainer&& x) {
  for (size_t i = 0; i < x.size(); ++i) {
    dest = save_varis(dest, x[i]);
  }
  return dest;
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
template <typename EigT, require_eigen_vt<is_var, EigT>* = nullptr>
inline vari** save_varis(vari** dest, EigT&& x) {
  for (int i = 0; i < x.size(); ++i) {
    dest[i] = x.coeffRef(i).vi_;
  }
  return dest + x.size();
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
template <typename Arith, require_st_arithmetic<Arith>* = nullptr>
inline vari** save_varis(vari** dest, Arith&& x) {
  return dest;
}

/**
 * End save_varis recursion and return pointer
 *
 * @param dest Pointer
 */
inline vari** save_varis(vari** dest) { return dest; }

}  // namespace math
}  // namespace stan

#endif
