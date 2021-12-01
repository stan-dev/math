#ifndef STAN_MATH_OPENCL_REV_DOT_PRODUCT_HPP
#define STAN_MATH_OPENCL_REV_DOT_PRODUCT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product.
 *
 * @tparam T1 type of the first vector
 * @tparam T2 type of the second vector
 *
 * @param[in] v1 First vector.
 * @param[in] v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error if sizes of v1 and v2 do not match.
 */
template <
    typename T1, typename T2, require_any_var_t<T1, T2>* = nullptr,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T1, T2>* = nullptr>
inline var dot_product(T1&& v1, T2&& v2) {
  arena_t<T1> v1_arena = std::forward<T1>(v1);
  arena_t<T2> v2_arena = std::forward<T2>(v2);

  return make_callback_var(dot_product(value_of(v1_arena), value_of(v2_arena)),
                           [v1_arena, v2_arena](vari& res) mutable {
                             adjoint_results(v1_arena, v2_arena)
                                 += expressions(res.adj() * value_of(v2_arena),
                                                res.adj() * value_of(v1_arena));
                           });
}

}  // namespace math
}  // namespace stan

#endif
#endif
