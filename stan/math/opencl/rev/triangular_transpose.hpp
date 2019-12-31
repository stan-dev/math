#ifndef STAN_MATH_OPENCL_REV_TRIANGULAR_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_REV_TRIANGULAR_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/triangular_transpose.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>

namespace stan {
namespace math {

/**
 * Copies a lower/upper triangular of a matrix to it's upper/lower.
 *
 * @tparam triangular_map Specifies if the copy is
 * lower-to-upper or upper-to-lower triangular. The value
 * must be of type TriangularMap
 *
 * @throw <code>std::invalid_argument</code> if the matrix is not square.
 *
 */
template <typename T>
template <TriangularMapCL triangular_map>
inline void matrix_cl<T, require_var_t<T>>::triangular_transpose() try {
  this->val().template triangular_transpose<triangular_map>();
  this->adj().template triangular_transpose<triangular_map>();
} catch (const cl::Error& e) {
  check_opencl_error("triangular_transpose", e);
}
}  // namespace math
}  // namespace stan

#endif
#endif
