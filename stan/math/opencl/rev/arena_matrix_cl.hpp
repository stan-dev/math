#ifndef STAN_MATH_OPENCL_REV_ARENA_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_REV_ARENA_MATRIX_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <utility>

namespace stan {
namespace math {
namespace internal {
template <typename T>
class arena_matrix_cl_impl : public chainable_alloc, public matrix_cl<T> {
 public:
  using Scalar = typename matrix_cl<T>::Scalar;
  using type = typename matrix_cl<T>::type;

  template <typename... Args>
  explicit arena_matrix_cl_impl(Args&&... args)
      : chainable_alloc(), matrix_cl<T>(std::forward<Args>(args)...) {}

  arena_matrix_cl_impl(const arena_matrix_cl_impl<T>&) = default;
  arena_matrix_cl_impl(arena_matrix_cl_impl<T>&) = default;
  arena_matrix_cl_impl(arena_matrix_cl_impl<T>&&) = default;
  arena_matrix_cl_impl<T>& operator=(const arena_matrix_cl_impl<T>&) = default;
  arena_matrix_cl_impl<T>& operator=(arena_matrix_cl_impl<T>&&) = default;
  using matrix_cl<T>::operator=;
};

}  // namespace internal

/**
 * A variant of `matrix_cl` that schedules its destructor to be called, so it
 * can be used on the AD stack.
 */
template <typename T>
class arena_matrix_cl : public matrix_cl_base {
 private:
  internal::arena_matrix_cl_impl<T>* impl_;
  template <typename>
  friend class matrix_cl;

 public:
  using Scalar = typename matrix_cl<T>::Scalar;
  using type = typename matrix_cl<T>::type;

  /**
   * General constructor forwards arguments to various `matrix_cl` constructors.
   * @tparam Args argument types
   * @param args arguments
   */
  template <typename... Args>
  explicit arena_matrix_cl(Args&&... args)
      : impl_(new internal::arena_matrix_cl_impl<T>(
            std::forward<Args>(args)...)) {}

  arena_matrix_cl(const arena_matrix_cl<T>&) = default;
  arena_matrix_cl(arena_matrix_cl<T>&) = default;
  arena_matrix_cl(arena_matrix_cl<T>&&) = default;
  arena_matrix_cl<T>& operator=(const arena_matrix_cl<T>&) = default;
  arena_matrix_cl<T>& operator=(arena_matrix_cl<T>&&) = default;

  /**
   * Constructor from a kernel generator expression.
   * @tparam Expr expression type
   * @param expression expression
   */
  // we need this as a separate overload, because the general constructor is
  // explicit
  template <typename Expr,
            require_all_kernel_expressions_and_none_scalar_t<Expr>* = nullptr>
  arena_matrix_cl(Expr&& expression)  // NOLINT(runtime/explicit)
      : impl_(new internal::arena_matrix_cl_impl<T>(
            std::forward<Expr>(expression))) {}

  /**
   * Implicit conversion operator to `matrix_cl`.
   * @return `matrix_cl` equivalent to `*this`. Returned matrix references the
   * same underlying buffer.
   */
  operator const matrix_cl<T>&() const {  // NOLINT(runtime/explicit)
    return *static_cast<const matrix_cl<T>*>(impl_);
  }
  operator matrix_cl<T>&() {  // NOLINT(runtime/explicit)
    return *static_cast<matrix_cl<T>*>(impl_);
  }

  /**
   * Evaluates `this`.
   * @return `matrix_cl` equivalent to `*this`
   */
  matrix_cl<T> eval() const& { return *impl_; }
  matrix_cl<T> eval() && { return std::move(*impl_); }

  // Wrapers to functions with explicit template parameters are implemented
  // without macros.
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros_strict_tri() {
    impl_->template zeros_strict_tri<matrix_view>();
  }

/**
 * Implements a wrapper for a non-const function in `matrix_cl`.
 * @param function_name name of the function to wrap
 */
#define ARENA_MATRIX_CL_FUNCTION_WRAPPER(function_name)       \
  template <typename... Args>                                 \
  inline decltype(auto) function_name(Args&&... args) {       \
    return impl_->function_name(std::forward<Args>(args)...); \
  }

  /**
   * Implements a wrapper for a const function in `matrix_cl`.
   * @param function_name name of the function to wrap
   */
#define ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(function_name) \
  template <typename... Args>                                 \
  inline decltype(auto) function_name(Args&&... args) const { \
    return impl_->function_name(std::forward<Args>(args)...); \
  }

  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(rows)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(cols)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(size)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(view)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(clear_write_events)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(clear_read_events)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(clear_read_write_events)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(write_events)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(read_events)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(read_write_events)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(add_read_event)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(add_write_event)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(add_read_write_event)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(wait_for_write_events)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(wait_for_read_events)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(wait_for_read_write_events)
  ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER(buffer)
  ARENA_MATRIX_CL_FUNCTION_WRAPPER(buffer)
  ARENA_MATRIX_CL_FUNCTION_WRAPPER(operator=)

#undef ARENA_MATRIX_CL_FUNCTION_WRAPPER
#undef ARENA_MATRIX_CL_CONST_FUNCTION_WRAPPER
};
template <typename T>
matrix_cl<T>::matrix_cl(const arena_matrix_cl<T>& A)
    // works like a move constructor, except it does not modify `a`
    : buffer_cl_(A.impl_->buffer_cl_),
      rows_(A.impl_->rows_),
      cols_(A.impl_->cols_),
      view_(A.impl_->view_),
      write_events_(A.impl_->write_events_),
      read_events_(A.impl_->read_events_) {}

template <typename T>
matrix_cl<T>& matrix_cl<T>::operator=(const arena_matrix_cl<T>& a) {
  // works like a move assignment operator, except it does not modify `a`
  view_ = a.impl_->view();
  rows_ = a.impl_->rows();
  cols_ = a.impl_->cols();
  this->wait_for_read_write_events();
  buffer_cl_ = a.impl_->buffer_cl_;
  write_events_ = a.impl_->write_events_;
  read_events_ = a.impl_->read_events_;
  return *this;
}

}  // namespace math
}  // namespace stan

#endif
#endif
