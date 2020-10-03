#ifndef STAN_MATH_REV_CORE_ACCUMULATE_ADJOINTS_HPP
#define STAN_MATH_REV_CORE_ACCUMULATE_ADJOINTS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>

#include <utility>
#include <vector>

namespace stan {
namespace math {

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
inline void accumulate_adjoints(double* dest, size_t pos, const var& x) {
  *(dest + pos) += x.adj();
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
template <typename VarVec, require_std_vector_vt<is_var, VarVec>* = nullptr>
inline void accumulate_adjoints(double* dest, size_t pos, VarVec&& x) {
  auto* new_dest = dest + pos;
  for (auto&& x_iter : x) {
    *new_dest += x_iter.adj();
    ++new_dest;
  }
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
template <typename VecContainer, require_std_vector_st<is_var, VecContainer>* = nullptr,
          require_std_vector_vt<is_container, VecContainer>* = nullptr>
inline void accumulate_adjoints(double* dest, size_t pos, VecContainer&& x) {
  size_t start_pos = pos;
  for (auto&& x_iter : x) {
    accumulate_adjoints(dest, start_pos, x_iter);
    start_pos += x_iter.size();
  }
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
template <typename EigT, require_eigen_vt<is_var, EigT>* = nullptr>
inline void accumulate_adjoints(double* dest, size_t pos, EigT&& x) {
  Eigen::Map<Eigen::MatrixXd>(dest + pos, x.rows(), x.cols()) += x.adj();
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
template <typename Arith, require_st_arithmetic<Arith>* = nullptr>
inline void accumulate_adjoints(double* dest, size_t pos, Arith&& x) {
}


}  // namespace math
}  // namespace stan

#endif
