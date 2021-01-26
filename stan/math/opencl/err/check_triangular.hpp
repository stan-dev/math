#ifndef STAN_MATH_OPENCL_ERR_CHECK_TRIANGULAR_HPP
#define STAN_MATH_OPENCL_ERR_CHECK_TRIANGULAR_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err/invalid_argument.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
namespace stan {
namespace math {

/** \ingroup error_checks_opencl
 * Check if the <code>matrix_cl</code> is either upper triangular or lower
 * triangular.
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param A <code>matrix_cl</code> to test
 *
 * @throw <code>std::invalid_argument</code> if the <code>matrix_cl</code>
 *    is not triangular
 */
template <typename T, require_matrix_cl_t<T>* = nullptr>
inline void check_triangular(const char* function, const char* name,
                             const T& A) {
  if (A.view() != matrix_cl_view::Lower && A.view() != matrix_cl_view::Upper) {
    invalid_argument("tri_inverse(OpenCL)", "A.view()",
                     static_cast<int>(A.view()), "is ",
                     ". Only triangular input matrices are supported!");
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
