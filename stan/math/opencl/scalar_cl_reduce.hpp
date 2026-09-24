#ifndef STAN_MATH_OPENCL_SCALAR_CL_REDUCE_HPP
#define STAN_MATH_OPENCL_SCALAR_CL_REDUCE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/kernels/scalar_reduce.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/scalar_cl.hpp>
#include <algorithm>
#include <type_traits>

namespace stan {
namespace math {
namespace opencl {
namespace internal {

/**
 * Work group size for the single work group reduction kernels: the largest
 * power of two no larger than 256 or the device's `LOCAL_SIZE_`.
 * @return work group size
 */
inline int scalar_reduce_local_size() {
  const int limit = std::min(256, opencl_context.base_opts().at("LOCAL_SIZE_"));
  int local_size = 1;
  while (local_size * 2 <= limit) {
    local_size *= 2;
  }
  return local_size;
}

/**
 * Reduces a non-empty kernel generator expression into a device scalar
 * without transferring data to the host. Large inputs are first reduced to
 * partial results with a 2D reduction; the partial results are then reduced by
 * a single work group.
 * @tparam Dst type of the device scalar or writable view receiving the result
 * @tparam T type of the expression
 * @tparam Partial callable computing a 2D partial reduction of an expression
 * @param[in,out] dst device scalar receiving the result
 * @param m expression to reduce
 * @param partial_f computes the 2D partial reduction
 * @param kernel single work group reduction kernel
 * @param accumulate whether to combine the result with `dst`
 * @param offset value added to the result
 */
template <typename Dst, typename T, typename Partial, typename Kernel>
inline void reduce_into(Dst&& dst, const T& m, Partial&& partial_f,
                        const Kernel& kernel, bool accumulate, double offset) {
  if constexpr (!is_matrix_cl<T>::value
                && std::decay_t<decltype(
                    as_operation_cl(m))>::Deriv::require_specific_local_size) {
    // an expression containing a reduction can not be nested in another
    // reduction, so it is evaluated first
    reduce_into(std::forward<Dst>(dst), matrix_cl<double>(m),
                std::forward<Partial>(partial_f), kernel, accumulate, offset);
    return;
  } else {
    const int local_size = scalar_reduce_local_size();
    const auto reduce = [&](const matrix_cl<double>& partials) {
      kernel(cl::NDRange(local_size), cl::NDRange(local_size), dst, partials,
             partials.size(), offset, static_cast<int>(accumulate));
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
      partials = partial_f(transpose(m));
    } else {
      partials = partial_f(m);
    }
    reduce(partials);
  }
}

/**
 * Sums a kernel generator expression into a device scalar without
 * transferring data to the host.
 * @tparam Dst type of the device scalar or writable view receiving the sum
 * @tparam T type of the expression
 * @param[in,out] dst device scalar receiving the sum
 * @param m expression to sum
 * @param accumulate whether to add the sum to `dst` instead of overwriting it
 * @param offset value added to the sum
 */
template <typename Dst, typename T, require_prim_scalar_cl_t<Dst>* = nullptr,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline void sum_into(Dst&& dst, const T& m, bool accumulate,
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
  reduce_into(
      std::forward<Dst>(dst), m, [](const auto& x) { return sum_2d(x); },
      opencl_kernels::scalar_sum, accumulate, offset);
}

/**
 * Computes the maximum of a kernel generator expression as a device scalar
 * without transferring data to the host. The maximum of no elements is
 * negative infinity.
 * @tparam T type of the expression
 * @param m expression
 * @return device scalar holding the maximum
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline ScalarCl<double> max_scalar_cl(const T& m) {
  if (m.size() == 0) {
    return ScalarCl<double>(scalar_<double>(NEGATIVE_INFTY));
  }
  ScalarCl<double> res(matrix_cl<double>(1, 1));
  reduce_into(
      res, m, [](const auto& x) { return max_2d(x); },
      opencl_kernels::scalar_max, false, 0.0);
  return res;
}

/**
 * Computes the product of a kernel generator expression as a device scalar
 * without transferring data to the host. The product of no elements is one.
 * @tparam T type of the expression
 * @param m expression
 * @return device scalar holding the product
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline ScalarCl<double> prod_scalar_cl(const T& m) {
  if (m.size() == 0) {
    return ScalarCl<double>(scalar_<double>(1.0));
  }
  ScalarCl<double> res(matrix_cl<double>(1, 1));
  reduce_into(
      res, m, [](const auto& x) { return prod_2d(x); },
      opencl_kernels::scalar_prod, false, 0.0);
  return res;
}

}  // namespace internal
}  // namespace opencl
}  // namespace math
}  // namespace stan

#endif
#endif
