#ifndef STAN_MATH_OPENCL_PRIM_PROD_HPP
#define STAN_MATH_OPENCL_PRIM_PROD_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/scalar_cl_reduce.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Calculates product of given kernel generator expression elements.
 *
 * Floating point products are computed on the device and returned as a device
 * scalar; no data is transferred to the host. Integer products are returned as
 * a host value.
 * @tparam T type of the expression
 * @param m expression to calculate product of
 * @return product of given expression, `opencl::ScalarCl<double>` for floating
 * point expressions
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline auto prod(const T& m) {
  if constexpr (std::is_floating_point<value_type_t<T>>::value) {
    return opencl::internal::prod_scalar_cl(m);
  } else {
    if constexpr (is_matrix_cl<T>::value) {
      if (m.size() < 1000) {
        // for small matrices running another kernel is not worth it
        return prod(from_matrix_cl(m));
      }
    }
    matrix_cl<value_type_t<T>> res;
    if (m.rows() <= 8) {
      // without transpose we would use just a few threads in a work group
      res = prod_2d(transpose(m));
    } else {
      res = prod_2d(m);
    }
    return prod(from_matrix_cl(res));
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
