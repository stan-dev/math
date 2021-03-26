#ifndef STAN_MATH_REV_CORE_DEEP_COPY_VARS_HPP
#define STAN_MATH_REV_CORE_DEEP_COPY_VARS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>

#include <utility>
#include <vector>

namespace stan {
namespace math {

/**
 * Forward arguments that do not contain vars. There
 *   is no copying to be done.
 *
 * @tparam Arith an arithmetic type.
 * @param arg For lvalue references this will be passed by reference.
 *  Otherwise it will be moved.
 */
template <typename Arith, typename = require_arithmetic_t<scalar_type_t<Arith>>>
inline Arith deep_copy_vars(Arith&& arg) {
  return std::forward<Arith>(arg);
}

/**
 * Copy the value of a var but reallocate a new vari
 *
 * @param arg A var
 * @return A new var
 */
inline auto deep_copy_vars(const var& arg) {
  return var(new vari(arg.val(), false));
}

/**
 * Copy the vars in arg but reallocate new varis for them
 *
 * @tparam VarVec A variant of std::vector<var>
 * @param arg A std::vector of vars
 * @return A new std::vector of vars
 */
template <typename VarVec, require_std_vector_vt<is_var, VarVec>* = nullptr>
inline auto deep_copy_vars(VarVec&& arg) {
  std::vector<var> copy_vec(arg.size());
  for (size_t i = 0; i < arg.size(); ++i) {
    copy_vec[i] = new vari(arg[i].val(), false);
  }
  return copy_vec;
}

/**
 * Copy the vars in arg but reallocate new varis for them
 *
 * @tparam VecContainer std::vector<T> where T is another type containing vars
 * @param arg A std::vector of containers containing vars
 * @return A new std::vector of containers containing vars
 */
template <typename VecContainer,
          require_std_vector_st<is_var, VecContainer>* = nullptr,
          require_std_vector_vt<is_container, VecContainer>* = nullptr>
inline auto deep_copy_vars(VecContainer&& arg) {
  std::vector<value_type_t<VecContainer>> copy_vec(arg.size());
  for (size_t i = 0; i < arg.size(); ++i) {
    copy_vec[i] = deep_copy_vars(arg[i]);
  }
  return copy_vec;
}

/**
 * Copy the vars in arg but reallocate new varis for them
 *
 * @tparam EigT An Eigen type with var value type
 * @param arg An Eigen container of vars
 * @return A new Eigen container of vars
 */
template <typename EigT, require_eigen_vt<is_var, EigT>* = nullptr>
inline auto deep_copy_vars(EigT&& arg) {
  return arg.unaryExpr([](auto&& x) { return var(new vari(x.val(), false)); })
      .eval();
}

}  // namespace math
}  // namespace stan

#endif
