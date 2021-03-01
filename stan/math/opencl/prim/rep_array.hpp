#ifndef STAN_MATH_OPENCL_PRIM_REP_ARRAY_HPP
#define STAN_MATH_OPENCL_PRIM_REP_ARRAY_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/rep_matrix.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Creates a matrix_cl representing an array
 * by replicating the input value.
 *
 * @tparam T type of the input
 * @param x the input value
 * @param n number of elements in the result array
 *
 * @return matrix_cl with replicated value from the input matrix
 *
 * @throw <code>domain_error</code> if the
 * requested dimensions are negative
 *
 */
template <
    typename T,
    require_any_t<
        is_matrix_cl<T>,
        math::conjunction<is_var<T>, is_matrix_cl<value_type_t<T>>>>* = nullptr>
auto rep_array(const scalar_type_t<T>& x, int n) {
  return rep_matrix<T>(x, n, 1);
}

}  // namespace math
}  // namespace stan

#endif
#endif
