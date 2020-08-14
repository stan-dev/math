#ifndef STAN_MATH_OPENCL_REV_ARENA_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_REV_ARENA_MATRIX_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * A variant of `matrix_cl` that schedules its destructor to be called, so it
 * can be used on the AD stack.
 */
template <typename T>
class arena_matrix_cl : public chainable_alloc, public matrix_cl<T> {
 public:
  using Scalar = typename matrix_cl<T>::Scalar;
  using type = typename matrix_cl<T>::type;

  /**
   * General constructoer forwards arguments to various `matrix_cl` constructors
   * and schedules the destructor to be called.
   * @tparam Args argument types
   * @param args arguments
   */
  template <typename... Args>
  explicit arena_matrix_cl(Args&&... args)
      : chainable_alloc(), matrix_cl<T>(std::forward<Args>(args)...) {}

  arena_matrix_cl(const arena_matrix_cl&) = default;
  arena_matrix_cl(arena_matrix_cl&) = default;
  arena_matrix_cl(arena_matrix_cl&&) = default;

  /**
   * Constructor from a kernel generator expression.
   * @tparam Expr expression type
   * @param expression expression
   */
  template <typename Expr,
            require_all_kernel_expressions_and_none_scalar_t<Expr>* = nullptr>
  arena_matrix_cl(const Expr& expression)
      : chainable_alloc(),
        matrix_cl<T>(expression) {}  // NOLINT(runtime/explicit)
};

}  // namespace math
}  // namespace stan

#endif
#endif
