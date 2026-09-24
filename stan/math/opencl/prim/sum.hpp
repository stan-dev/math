#ifndef STAN_MATH_OPENCL_PRIM_SUM_HPP
#define STAN_MATH_OPENCL_PRIM_SUM_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/kernels/scalar_sum.hpp>
#include <stan/math/opencl/scalar_cl.hpp>
#include <algorithm>

namespace stan {
namespace math {

namespace opencl {
namespace internal {

/**
 * Sums a kernel generator expression into a device scalar without
 * transferring data to the host. Large inputs are first reduced to partial
 * sums with `sum_2d()`; the partial sums are then reduced by a single work
 * group.
 * @tparam T type of the expression
 * @param[in,out] dst device scalar receiving the sum
 * @param m expression to sum
 * @param accumulate whether to add the sum to `dst` instead of overwriting it
 * @param offset value added to the sum
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline void sum_into(ScalarCl<double>& dst, const T& m, bool accumulate,
                     double offset = 0.0) {
  if (m.size() == 0) {
    if (accumulate) {
      if (offset != 0.0) {
        dst += offset;
      }
    } else {
      dst = scalar_<double>(offset);
    }
    return;
  }
  // The reduction kernel requires a power of two work group size.
  int local_size = std::min(256, opencl_context.base_opts().at("LOCAL_SIZE_"));
  int pow2 = 1;
  while (pow2 * 2 <= local_size) {
    pow2 *= 2;
  }
  local_size = pow2;
  const auto reduce = [&](const matrix_cl<double>& partials) {
    opencl_kernels::scalar_sum(cl::NDRange(local_size), cl::NDRange(local_size),
                               dst, partials, partials.size(), offset,
                               static_cast<int>(accumulate));
  };
  if constexpr (is_matrix_cl<T>::value) {
    // A buffer with no implied zeros can be reduced directly.
    if (m.view() == matrix_cl_view::Entire && m.size() <= 64 * local_size) {
      reduce(m);
      return;
    }
  }
  matrix_cl<double> partials;
  if (m.rows() <= 8) {
    // without transpose we would use just a few threads in a work group
    partials = sum_2d(transpose(m));
  } else {
    partials = sum_2d(m);
  }
  reduce(partials);
}

}  // namespace internal
}  // namespace opencl

/**
 * Calculates sum of given kernel generator expression.
 * @tparam T type of the expression
 * @param m expression to sum
 * @return sum of given expression
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline value_type_t<T> sum(const T& m) {
  if constexpr (is_matrix_cl<T>::value) {
    if (m.size() < 1000) {
      // for small matrices running another kernel is not worth it
      return sum(from_matrix_cl(m));
    }
  }
  matrix_cl<value_type_t<T>> res;
  if (m.rows() <= 8) {
    // without transpose we would use just a few threads in a work group
    res = sum_2d(transpose(m));
  } else {
    res = sum_2d(m);
  }
  return sum(from_matrix_cl(res));
}

}  // namespace math
}  // namespace stan

#endif
#endif
