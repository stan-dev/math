#ifndef STAN_MATH_REV_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_REV_FUN_VALUE_OF_REC_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * @param v Variable.
 * @return Value of variable.
 */
template <typename Var, require_var_t<Var>* = nullptr>
inline auto value_of_rec(Var&& v) {
  return v.vi_->val_;
}

template <typename Vec, require_std_vector_vt<is_var, Vec>* = nullptr>
inline auto value_of_rec(Vec&& x) {
  return value_of(std::forward<Vec>(x));
}

template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
inline auto value_of_rec(EigMat&& M) {
  return M.val().eval();
}

}  // namespace math
}  // namespace stan
#endif
