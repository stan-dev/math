#ifndef STAN_MATH_REV_FUN_VALUE_OF_HPP
#define STAN_MATH_REV_FUN_VALUE_OF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * <p>This function is used internally by autodiff functions along
 * with <code>value_of(T x)</code> to extract the
 * <code>double</code> value of either a scalar or an autodiff
 * variable.  This function will be called when the argument is a
 * <code>var</code> even if the function is not
 * referred to by namespace because of argument-dependent lookup.
 *
 * @param v Variable.
 * @return Value of variable.
 */
template <typename T>
inline auto& value_of(const var_value<T>& v) {
  return v.vi_->val_;
}

template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
inline auto value_of(EigMat&& v) {
  return make_holder([](auto&& a) {
    return std::forward<decltype(a)>(a).val();
  }, std::forward<EigMat>(v));
}

}  // namespace math
}  // namespace stan
#endif
