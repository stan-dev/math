#ifndef STAN_MATH_REV_CORE_DEEP_COPY_VARS_HPP
#define STAN_MATH_REV_CORE_DEEP_COPY_VARS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor.hpp>
#include <stan/math/rev/fun.hpp>

#include <iostream>
#include <iterator>
#include <tuple>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

/**
 * Specialization of deep copy that returns references for arithmetic types.
 * @tparam Arith an arithmetic type.
 * @param arg For lvalue references this will be passed by reference.
 *  Otherwise it will be moved.
 */
template <typename Arith, typename = require_arithmetic_t<scalar_type_t<Arith>>>
inline decltype(auto) deep_copy_vars(Arith&& arg) {
  return std::forward<Arith>(arg);
}

/**
 * Specialization to copy a single var
 * @param arg a var.
 */
inline auto deep_copy_vars(const var& arg) {
  return var(new vari(arg.val(), false));
}

/**
 * Specialization to copy a standard vector of vars.
 * @tparam VarVec A standard vector holding vars.
 * @param arg A standard vector holding vars.
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
 * Specialization to copy a standard vector holding containers.
 * @tparam VarVec A standard vector holding containers of vars.
 * @param arg A standard vector holding containers of vars.
 */
template <typename VarVec, require_std_vector_st<is_var, VarVec>* = nullptr,
          require_std_vector_vt<is_container, VarVec>* = nullptr>
inline auto deep_copy_vars(VarVec&& arg) {
  std::vector<value_type_t<VarVec>> copy_vec(arg.size());
  for (size_t i = 0; i < arg.size(); ++i) {
    copy_vec[i] = deep_copy_vars(arg[i]);
  }
  return copy_vec;
}

/**
 * Specialization to copy an `Eigen` containers of vars.
 * @tparam EigT A type derived from `EigenBase` containing vars.
 * @param arg A container hollding var types.
 */
template <typename EigT, require_eigen_vt<is_var, EigT>* = nullptr>
inline auto deep_copy_vars(EigT&& arg) {
  return arg.unaryExpr([](auto&& x) { return var(new vari(x.val(), false)); })
      .eval();
}

}  // namespace math
}  // namespace stan

#endif
