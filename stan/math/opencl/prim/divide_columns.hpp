#ifndef STAN_MATH_OPENCL_PRIM_DIVIDE_COLUMNS_HPP
#define STAN_MATH_OPENCL_PRIM_DIVIDE_COLUMNS_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/divide_columns.hpp>
#include <stan/math/opencl/err.hpp>
#include <CL/opencl.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Divides each column of a matrix by a vector
 *
 * @tparam T1 Type of first matrix
 * @tparam T2 Type of second matrix
 * @param A first matrix
 * @param B Vector of elements to divide each column of \c A by element-wise.
 *
 * @return element-wise division of \c A by \c B.
 *
 * @throw <code>std::invalid_argument</code> if the
 * input matrices do not have matching dimensions
 *
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline void divide_columns(matrix_cl<T1>& A, const matrix_cl<T2>& B) {
  if (A.size() == 0 || B.size() == 0) {
    return;
  }
  check_size_match("divide_columns", "A mod B", A.size() % B.size(), "B mod",
                   0);
  check_vector("divide_columns", "B", B);
  try {
    opencl_kernels::divide_columns_vec(cl::NDRange(A.size()), A, B, B.size());
  } catch (const cl::Error& e) {
    check_opencl_error("divide_columns", e);
  }
}

/** \ingroup opencl
 * Divides each column of a matrix by a scalar
 *
 * @tparam T1 Type of first matrix
 * @tparam T2 Type of the divisor
 * @param A Matrix to divide
 * @param divisor scalar to divide by
 *
 * @return element-wise division of \c A by \c divisor.
 *
 */
template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline void divide_columns(matrix_cl<T1>& A, const T2& divisor) {
  A = elt_divide(A, divisor);
}

}  // namespace math
}  // namespace stan

#endif
#endif
