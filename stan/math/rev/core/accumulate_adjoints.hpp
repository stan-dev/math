#ifndef STAN_MATH_REV_CORE_ACCUMULATE_ADJOINTS_HPP
#define STAN_MATH_REV_CORE_ACCUMULATE_ADJOINTS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>

#include <utility>
#include <vector>

namespace stan {
namespace math {

template <typename... Pargs>
inline double* accumulate_adjoints(double* dest, const var& x, Pargs&&... args);

template <typename VarVec, require_std_vector_vt<is_var, VarVec>* = nullptr,
          typename... Pargs>
inline double* accumulate_adjoints(double* dest, VarVec&& x, Pargs&&... args);

template <typename VecContainer,
          require_std_vector_st<is_var, VecContainer>* = nullptr,
          require_std_vector_vt<is_container, VecContainer>* = nullptr,
          typename... Pargs>
inline double* accumulate_adjoints(double* dest, VecContainer&& x,
                                   Pargs&&... args);

template <typename EigT, require_eigen_vt<is_var, EigT>* = nullptr,
          typename... Pargs>
inline double* accumulate_adjoints(double* dest, EigT&& x, Pargs&&... args);

template <typename Arith, require_st_arithmetic<Arith>* = nullptr,
          typename... Pargs>
inline double* accumulate_adjoints(double* dest, Arith&& x, Pargs&&... args);

inline double* accumulate_adjoints(double* dest);

/**
 * Accumulate adjoints from x into storage pointed to by dest,
 *   increment the adjoint storage pointer,
 *   recursively accumulate the adjoints of the rest of the arguments,
 *   and return final position of storage pointer.
 *
 * @tparam Pargs Types of remaining arguments
 * @param dest Pointer to where adjoints are to be accumulated
 * @param x A var
 * @param args Further args to accumulate over
 * @return Final position of adjoint storage pointer
 */
template <typename... Pargs>
inline double* accumulate_adjoints(double* dest, const var& x,
                                   Pargs&&... args) {
  *dest += x.adj();
  return accumulate_adjoints(dest + 1, std::forward<Pargs>(args)...);
}

/**
 * Accumulate adjoints from std::vector x into storage pointed to by dest,
 *   increment the adjoint storage pointer,
 *   recursively accumulate the adjoints of the rest of the arguments,
 *   and return final position of storage pointer.
 *
 * @tparam Pargs Types of remaining arguments
 * @param dest Pointer to where adjoints are to be accumulated
 * @param x A std::vector of vars
 * @param args Further args to accumulate over
 * @return Final position of adjoint storage pointer
 */
template <typename VarVec, require_std_vector_vt<is_var, VarVec>*,
          typename... Pargs>
inline double* accumulate_adjoints(double* dest, VarVec&& x, Pargs&&... args) {
  for (auto&& x_iter : x) {
    *dest += x_iter.adj();
    ++dest;
  }
  return accumulate_adjoints(dest, std::forward<Pargs>(args)...);
}

/**
 * Accumulate adjoints from x (a std::vector of containers containing vars)
 *   into storage pointed to by dest,
 *   increment the adjoint storage pointer,
 *   recursively accumulate the adjoints of the rest of the arguments,
 *   and return final position of storage pointer.
 *
 * @tparam VecContainer the type of a standard container holding var
 * containers.
 * @tparam Pargs Types of remaining arguments
 * @param dest Pointer to where adjoints are to be accumulated
 * @param x A std::vector of containers holding vars
 * @param args Further args to accumulate over
 * @return Final position of adjoint storage pointer
 */
template <typename VecContainer, require_std_vector_st<is_var, VecContainer>*,
          require_std_vector_vt<is_container, VecContainer>*, typename... Pargs>
inline double* accumulate_adjoints(double* dest, VecContainer&& x,
                                   Pargs&&... args) {
  for (auto&& x_iter : x) {
    dest = accumulate_adjoints(dest, x_iter);
  }
  return accumulate_adjoints(dest, std::forward<Pargs>(args)...);
}

/**
 * Accumulate adjoints from x (an Eigen type containing vars)
 *   into storage pointed to by dest,
 *   increment the adjoint storage pointer,
 *   recursively accumulate the adjoints of the rest of the arguments,
 *   and return final position of storage pointer.
 *
 * @tparam EigT Type derived from `EigenBase` containing vars.
 * @tparam Pargs Types of remaining arguments
 * @param dest Pointer to where adjoints are to be accumulated
 * @param x An eigen type holding vars to accumulate over
 * @param args Further args to accumulate over
 * @return Final position of adjoint storage pointer
 */
template <typename EigT, require_eigen_vt<is_var, EigT>*, typename... Pargs>
inline double* accumulate_adjoints(double* dest, EigT&& x, Pargs&&... args) {
  Eigen::Map<Eigen::MatrixXd>(dest, x.rows(), x.cols()) += x.adj();
  return accumulate_adjoints(dest + x.size(), std::forward<Pargs>(args)...);
}

/**
 * Ignore arithmetic types.
 *
 * Recursively accumulate the adjoints of the rest of the arguments
 *   and return final position of adjoint storage pointer.
 *
 * @tparam Arith A type satisfying `std::is_arithmetic`.
 * @tparam Pargs Types of remaining arguments
 * @param dest Pointer to where adjoints are to be accumulated
 * @param x An object that is either arithmetic or a container of Arithmetic
 *  types
 * @param args Further args to accumulate over
 * @return Final position of adjoint storage pointer
 */
template <typename Arith, require_st_arithmetic<Arith>*, typename... Pargs>
inline double* accumulate_adjoints(double* dest, Arith&& x, Pargs&&... args) {
  return accumulate_adjoints(dest, std::forward<Pargs>(args)...);
}

/**
 * End accumulate_adjoints recursion and return pointer
 *
 * @param dest Pointer
 */
inline double* accumulate_adjoints(double* dest) { return dest; }

}  // namespace math
}  // namespace stan

#endif
