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
template <typename T, require_var_vt<is_matrix_cl, T>* = nullptr>
inline var dot_self(T&& v) {
  arena_t<T> v_arena = std::forward<T>(v);

  return make_callback_var(dot_self(v.val()), [v_arena](vari& res) mutable {
    v_arena.adj() = v_arena.adj() + 2.0 * res.adj() * v_arena.val();
  });
}

}  // namespace math
}  // namespace stan

#endif
#endif
