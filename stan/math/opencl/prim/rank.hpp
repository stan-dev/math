#ifndef STAN_MATH_OPENCL_PRIM_RANK_HPP
#define STAN_MATH_OPENCL_PRIM_RANK_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return the number of components of v less than v[s].
 *
 * @tparam C container type
 * @param[in] v input vector
 * @param[in] s position in vector
 * @return number of components of v less than v[s].
 * @throw std::out_of_range if s is out of range.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
int rank(const T& v, int s) {
  check_vector("rank (OpenCL)", "v", v);
  check_range("rank (OpenCL)", "v", v.size(), s);
  auto v_col = as_column_vector_or_scalar(v);
  return sum((v_col < indexing(v_col, s - 1, col_index())) + 0);
}

}  // namespace math
}  // namespace stan

#endif
#endif
