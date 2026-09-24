#ifndef STAN_MATH_OPENCL_PRIM_PARTIALS_PROPAGATOR_HPP
#define STAN_MATH_OPENCL_PRIM_PARTIALS_PROPAGATOR_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/scalar_cl.hpp>
#include <stan/math/opencl/scalar_cl_reduce.hpp>
#include <tuple>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
namespace opencl {

/** \ingroup opencl
 * The device scalar type returned by OpenCL functions of the given argument
 * types: `ScalarCl<var>` if any argument is autodiff, `ScalarCl<double>`
 * otherwise.
 */
template <typename... Ts>
using scalar_cl_return_t = ScalarCl<return_type_t<Ts...>>;

namespace internal {

/**
 * The partial derivative of a device result with respect to a scalar operand,
 * stored on the device.
 *
 * Assigning a matrix or a kernel generator expression sums it on the device,
 * which is the derivative with respect to a scalar that was broadcast over the
 * expression. Assigning a device scalar copies it. Host values must be moved
 * to the device explicitly.
 */
class scalar_cl_partial {
 public:
  ScalarCl<double> value_;

  /**
   * Assigns a device scalar or the sum of a kernel generator expression.
   * @tparam T type of the device scalar or expression
   * @param x device scalar or expression
   * @return this partial
   */
  template <typename T, require_all_kernel_expressions_t<T>* = nullptr,
            require_not_arithmetic_t<T>* = nullptr>
  scalar_cl_partial& operator=(T&& x) {
    if constexpr (is_prim_scalar_cl<T>::value) {
      value_ = ScalarCl<double>(x);
    } else {
      auto&& x_op = as_operation_cl(std::forward<T>(x));
      if (x_op.rows() == -1 && x_op.cols() == -1) {
        value_ = x_op;
      } else {
        sum_into(value_, x_op, false);
      }
    }
    return *this;
  }

  /**
   * Adds a device scalar or the sum of a kernel generator expression.
   * @tparam T type of the device scalar or expression
   * @param x device scalar or expression
   * @return this partial
   */
  template <typename T, require_all_kernel_expressions_t<T>* = nullptr,
            require_not_arithmetic_t<T>* = nullptr>
  scalar_cl_partial& operator+=(T&& x) {
    if constexpr (is_prim_scalar_cl<T>::value) {
      value_ += x;
    } else {
      auto&& x_op = as_operation_cl(std::forward<T>(x));
      if (x_op.rows() == -1 && x_op.cols() == -1) {
        value_ += x_op;
      } else {
        sum_into(value_, x_op, true);
      }
    }
    return *this;
  }

  /**
   * Multivariate functions index scalar partials; there is one element.
   * @return this partial
   */
  inline scalar_cl_partial& operator[](int) { return *this; }
};

}  // namespace internal
}  // namespace opencl

namespace internal {

/**
 * Partials propagator for OpenCL functions none of whose operands is
 * autodiff. Nothing is propagated and the value is returned as it is.
 */
template <typename ReturnType, typename... Ops>
class partials_propagator<ReturnType, require_prim_scalar_cl_t<ReturnType>,
                          Ops...> {
 public:
  template <typename... Types>
  explicit partials_propagator(Types&&... /* ops */) noexcept {}

  /**
   * @param value the value of the function
   * @return `value`
   */
  inline static opencl::ScalarCl<double> build(
      opencl::ScalarCl<double>&& value) noexcept {
    return std::move(value);
  }
};

}  // namespace internal

namespace opencl {

/** \ingroup opencl
 * Constructs a partials propagator for an OpenCL function. Its result is a
 * device scalar and derivatives stay on the device, whatever the operands are:
 * host scalars, CPU vars, device scalars or OpenCL matrices.
 *
 * @tparam Ops types of the operands
 * @param ops operands
 * @return partials propagator whose `build()` takes the device value of the
 * function
 */
template <typename... Ops>
inline auto make_partials_propagator(Ops&&... ops) {
  using return_type = scalar_cl_return_t<Ops...>;
  return math::internal::partials_propagator<
      return_type, void, plain_type_t<std::decay_t<Ops>>...>(
      std::forward<Ops>(ops)...);
}

}  // namespace opencl

}  // namespace math
}  // namespace stan

#endif
#endif
