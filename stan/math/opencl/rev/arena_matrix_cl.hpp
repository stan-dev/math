#ifndef STAN_MATH_OPENCL_REV_ARENA_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_REV_ARENA_MATRIX_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/core/needs_destructor.hpp>
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
class arena_matrix_cl : public needs_destructor, public matrix_cl<T> {
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
      : needs_destructor(), matrix_cl<T>(std::forward<Args>(args)...) {}

  arena_matrix_cl(const arena_matrix_cl& other) = default;
//      : matrix_cl<T>(other.buffer(), other.rows(), other.cols(), other.view()) {
//    read_events_ = other.read_events();
//    write_events_ = other.write_events();
//  }
  arena_matrix_cl(arena_matrix_cl&&) = default;

  /**
   * Constructor from a kernel generator expression.
   * @tparam Expr expression type
   * @param expression expression
   */
  template <typename Expr,
            require_all_kernel_expressions_and_none_scalar_t<Expr>* = nullptr>
  arena_matrix_cl(const Expr& expression)  // NOLINT(runtime/explicit)
      : needs_destructor(), matrix_cl<T>(expression) {}

//  void destroy() { ~matrix_cl<T>(); }
//  ~arena_matrix_cl() {
//    read_events_.~vector();
//    write_events_.~vector();
//  }
};

}  // namespace math
}  // namespace stan

#endif
#endif
