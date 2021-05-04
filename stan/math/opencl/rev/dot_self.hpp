#ifndef STAN_MATH_OPENCL_REV_DOT_SELF_HPP
#define STAN_MATH_OPENCL_REV_DOT_SELF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/dot_self.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of a vector of var with itself.
 *
 * @tparam T type of the vector
 * @param[in] v Vector.
 * @return Dot product of the vector with itself.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var dot_self(const var_value<T>& v) {
  return make_callback_var(dot_self(v.val()), [v](vari& res) mutable {
    v.adj() += 2.0 * res.adj() * v.val();
  });
}

}  // namespace math
}  // namespace stan

#endif
#endif
