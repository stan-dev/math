#ifndef STAN_MATH_OPENCL_SCALAR_CL_HPP
#define STAN_MATH_OPENCL_SCALAR_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <CL/opencl.hpp>
#include <tbb/concurrent_vector.h>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
namespace opencl {
namespace internal {

/**
 * Buffer and event access shared by device scalars and device scalar views.
 * Everything is forwarded to the backing 1x1 `matrix_cl` returned by
 * `Derived::matrix()`, so views and owners share event state.
 * @tparam Derived the device scalar type
 */
template <typename Derived>
class scalar_cl_access {
 public:
  /**
   * @return the OpenCL buffer holding the value
   */
  inline const cl::Buffer& buffer() const noexcept {
    return derived().matrix().buffer();
  }
  /**
   * @return events of all operations writing to the value
   */
  inline const tbb::concurrent_vector<cl::Event>& write_events() const {
    return derived().matrix().write_events();
  }
  /**
   * @return events of all operations reading the value
   */
  inline const tbb::concurrent_vector<cl::Event>& read_events() const {
    return derived().matrix().read_events();
  }
  /**
   * @return events of all operations reading or writing the value
   */
  inline tbb::concurrent_vector<cl::Event> read_write_events() const {
    return derived().matrix().read_write_events();
  }
  /**
   * Adds an event of an operation reading the value.
   * @param new_event event to add
   */
  inline void add_read_event(cl::Event new_event) const {
    derived().matrix().add_read_event(std::move(new_event));
  }
  /**
   * Adds an event of an operation writing the value.
   * @param new_event event to add
   */
  inline void add_write_event(cl::Event new_event) const {
    derived().matrix().add_write_event(std::move(new_event));
  }
  /**
   * Adds an event of an operation reading and writing the value.
   * @param new_event event to add
   */
  inline void add_read_write_event(cl::Event new_event) const {
    derived().matrix().add_read_write_event(std::move(new_event));
  }

 protected:
  inline const Derived& derived() const noexcept {
    return static_cast<const Derived&>(*this);
  }
  inline Derived& derived() noexcept { return static_cast<Derived&>(*this); }
};

/**
 * Assignment and compound assignment shared by writable device scalars. The
 * right hand side is evaluated on the device with a single thread.
 * @tparam Derived the device scalar type
 */
template <typename Derived>
class scalar_cl_assign : public scalar_cl_access<Derived> {
 public:
  /**
   * Adds a scalar or 1x1 kernel generator expression to this device scalar.
   * @tparam T type of the argument
   * @param b value to add
   * @return this device scalar
   */
  template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
  Derived& operator+=(T&& b) {
    assign(as_operation_cl(this->derived())
           + as_operation_cl(std::forward<T>(b)));
    return this->derived();
  }

  /**
   * Subtracts a scalar or 1x1 kernel generator expression from this device
   * scalar.
   * @tparam T type of the argument
   * @param b value to subtract
   * @return this device scalar
   */
  template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
  Derived& operator-=(T&& b) {
    assign(as_operation_cl(this->derived())
           - as_operation_cl(std::forward<T>(b)));
    return this->derived();
  }

  /**
   * Multiplies this device scalar by a scalar or 1x1 kernel generator
   * expression.
   * @tparam T type of the argument
   * @param b value to multiply by
   * @return this device scalar
   */
  template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
  Derived& operator*=(T&& b) {
    assign(elt_multiply(as_operation_cl(this->derived()),
                        as_operation_cl(std::forward<T>(b))));
    return this->derived();
  }

  /**
   * Divides this device scalar by a scalar or 1x1 kernel generator expression.
   * @tparam T type of the argument
   * @param b value to divide by
   * @return this device scalar
   */
  template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
  Derived& operator/=(T&& b) {
    assign(elt_divide(as_operation_cl(this->derived()),
                      as_operation_cl(std::forward<T>(b))));
    return this->derived();
  }

 protected:
  /**
   * Evaluates an expression into the backing buffer with a single thread.
   * @tparam Expr type of the expression
   * @param expression expression to evaluate
   * @throw std::invalid_argument if the expression is not 1x1
   */
  template <typename Expr>
  inline void assign(Expr&& expression) {
    this->derived().matrix() = scalar_result_<as_operation_cl_t<Expr>>(
        as_operation_cl(std::forward<Expr>(expression)));
  }
};

/**
 * Checks that a `matrix_cl` can back a device scalar.
 * @param m matrix to check
 * @throw std::invalid_argument if the matrix is not 1x1
 */
inline void check_scalar_cl_matrix(const matrix_cl<double>& m) {
  check_size_match("ScalarCl", "rows", m.rows(), "", 1);
  check_size_match("ScalarCl", "columns", m.cols(), "", 1);
}

}  // namespace internal

/** \ingroup opencl
 * A double that lives on the OpenCL device.
 *
 * The value is stored in a 1x1 `matrix_cl<double>` so it can be read and
 * written by kernels without being transferred to the host. There are no
 * implicit conversions to or from host values: use the explicit
 * constructor to upload a value and `opencl::to_host()` to read it back.
 *
 * Copies duplicate the value on the device. Moves transfer ownership of the
 * device buffer.
 */
template <>
class ScalarCl<double> : public internal::scalar_cl_assign<ScalarCl<double>> {
 private:
  matrix_cl<double> buf_;

 public:
  using value_type = double;

  /**
   * Constructs a device scalar holding zero. No host data is transferred.
   */
  ScalarCl() : buf_(1, 1) { buf_.setZero(); }

  /**
   * Uploads a host value to the device.
   *
   * The value is passed to `matrix_cl` as an rvalue so the write blocks until
   * it completes; an lvalue would enqueue a non-blocking write reading from
   * host memory that may no longer exist when the write executes.
   * @param value value to upload
   */
  explicit ScalarCl(double value)
      : buf_(std::move(value), matrix_cl_view::Entire) {}

  /**
   * Takes ownership of a 1x1 `matrix_cl`.
   * @param m matrix holding the value
   * @throw std::invalid_argument if the matrix is not 1x1
   */
  explicit ScalarCl(matrix_cl<double>&& m) : buf_(std::move(m)) {
    internal::check_scalar_cl_matrix(buf_);
  }

  /**
   * Copies the value of a device scalar view on the device.
   * @tparam T type of the view
   * @param other device scalar view
   */
  template <typename T, require_prim_scalar_cl_t<T>* = nullptr,
            require_not_same_t<std::decay_t<T>, ScalarCl<double>>* = nullptr>
  explicit ScalarCl(const T& other) : buf_(other.matrix()) {}

  /**
   * Evaluates a kernel generator expression into a device scalar. The
   * expression must be 1x1 or consist only of scalars.
   * @tparam Expr type of the expression
   * @param expression expression to evaluate
   * @throw std::invalid_argument if the expression is not 1x1
   */
  template <typename Expr,
            require_all_kernel_expressions_and_none_scalar_t<Expr>* = nullptr>
  explicit ScalarCl(Expr&& expression) : buf_(1, 1) {
    this->assign(std::forward<Expr>(expression));
  }

  ScalarCl(const ScalarCl& other) = default;
  ScalarCl(ScalarCl&& other) = default;
  ScalarCl& operator=(const ScalarCl& other) = default;
  ScalarCl& operator=(ScalarCl&& other) = default;

  /**
   * Evaluates a kernel generator expression into this device scalar. The
   * expression must be 1x1 or consist only of scalars.
   * @tparam Expr type of the expression
   * @param expression expression to evaluate
   * @return this device scalar
   * @throw std::invalid_argument if the expression is not 1x1
   */
  template <typename Expr,
            require_all_kernel_expressions_and_none_scalar_t<Expr>* = nullptr>
  ScalarCl& operator=(Expr&& expression) {
    this->assign(std::forward<Expr>(expression));
    return *this;
  }

  /**
   * @return the backing 1x1 `matrix_cl`
   */
  inline const matrix_cl<double>& matrix() const noexcept { return buf_; }
  /**
   * @return the backing 1x1 `matrix_cl`
   */
  inline matrix_cl<double>& matrix() noexcept { return buf_; }
};

/** \ingroup opencl
 * A writable view of a device scalar, such as the adjoint of a
 * `ScalarCl<var>`. It refers to a 1x1 `matrix_cl` owned elsewhere and shares
 * its buffer and events. Assigning to the view writes the referenced value.
 */
template <>
class ScalarCl<double&> : public internal::scalar_cl_assign<ScalarCl<double&>> {
 private:
  matrix_cl<double>* m_;

 public:
  using value_type = double;

  /**
   * @param m 1x1 matrix holding the value
   * @throw std::invalid_argument if the matrix is not 1x1
   */
  explicit ScalarCl(matrix_cl<double>& m) : m_(&m) {
    internal::check_scalar_cl_matrix(m);
  }

  ScalarCl(const ScalarCl& other) = default;

  /**
   * Assigns the value of another view to the value referenced by this one.
   * @param other view to copy the value from
   * @return this view
   */
  ScalarCl& operator=(const ScalarCl& other) {
    this->assign(other);
    return *this;
  }

  /**
   * Evaluates a device scalar or a kernel generator expression into the
   * referenced value. The expression must be 1x1 or consist only of scalars.
   * @tparam Expr type of the expression
   * @param expression expression to evaluate
   * @return this view
   * @throw std::invalid_argument if the expression is not 1x1
   */
  template <typename Expr, require_all_kernel_expressions_t<Expr>* = nullptr>
  ScalarCl& operator=(Expr&& expression) {
    this->assign(std::forward<Expr>(expression));
    return *this;
  }

  /**
   * @return the referenced 1x1 `matrix_cl`
   */
  inline matrix_cl<double>& matrix() const noexcept { return *m_; }
};

/** \ingroup opencl
 * A read-only view of a device scalar, such as the value of a
 * `ScalarCl<var>`. It refers to a 1x1 `matrix_cl` owned elsewhere and shares
 * its buffer and events.
 */
template <>
class ScalarCl<const double&>
    : public internal::scalar_cl_access<ScalarCl<const double&>> {
 private:
  const matrix_cl<double>* m_;

 public:
  using value_type = double;

  /**
   * @param m 1x1 matrix holding the value
   * @throw std::invalid_argument if the matrix is not 1x1
   */
  explicit ScalarCl(const matrix_cl<double>& m) : m_(&m) {
    internal::check_scalar_cl_matrix(m);
  }

  /**
   * Views the value of a device scalar or of a writable view.
   * @tparam T type of the device scalar
   * @param x device scalar
   */
  template <
      typename T, require_prim_scalar_cl_t<T>* = nullptr,
      require_not_same_t<std::decay_t<T>, ScalarCl<const double&>>* = nullptr>
  ScalarCl(const T& x)  // NOLINT(runtime/explicit)
      : m_(&x.matrix()) {}

  ScalarCl(const ScalarCl& other) = default;
  ScalarCl& operator=(const ScalarCl& other) = delete;

  /**
   * @return the referenced 1x1 `matrix_cl`
   */
  inline const matrix_cl<double>& matrix() const noexcept { return *m_; }
};

namespace internal {
/**
 * Prepares a value for a handwritten kernel that takes a buffer. A device
 * scalar passes its own 1x1 buffer, with no transfer; host values are copied
 * to the device.
 * @tparam T type of the value
 * @param x value
 * @return the backing `matrix_cl` of a device scalar, `to_matrix_cl(x)`
 * otherwise
 */
template <typename T>
inline decltype(auto) as_kernel_buffer(const T& x) {
  if constexpr (is_prim_scalar_cl<T>::value) {
    return x.matrix();
  } else {
    return to_matrix_cl(x);
  }
}

/**
 * Checks if a type is a scalar on the host or on the device, as opposed to a
 * container of values.
 */
template <typename T>
struct is_host_or_device_scalar
    : math::disjunction<is_stan_scalar<T>, is_scalar_cl<T>> {};

/**
 * Prepares a value for use in a kernel generator expression that also
 * involves matrices. Device scalars become a `scalar_buf_` operation, so
 * comparisons, checks and arithmetic on them are fused into the same kernel as
 * the matrix operations instead of being evaluated on their own. Other values
 * are returned unchanged.
 * @tparam T type of the value
 * @param x value
 * @return kernel generator operation for a device scalar, `x` otherwise
 */
template <typename T>
inline decltype(auto) as_operand(T&& x) {
  if constexpr (is_prim_scalar_cl<T>::value) {
    return as_operation_cl(std::forward<T>(x));
  } else {
    return std::forward<T>(x);
  }
}
}  // namespace internal

/** \ingroup opencl
 * A device scalar is already a scalar, so it is returned as it is.
 * @tparam T type of the device scalar
 * @param x device scalar
 * @return `x`
 */
template <typename T, require_scalar_cl_t<T>* = nullptr>
inline const T& as_column_vector_or_scalar(const T& x) {
  return x;
}

/** \ingroup opencl
 * Copies a device scalar to the host. Blocks until all writes to the scalar
 * have finished and the value has been read.
 * @tparam T type of the device scalar
 * @param x device scalar
 * @return host copy of the value
 */
template <typename T, require_prim_scalar_cl_t<T>* = nullptr>
inline double to_host(const T& x) {
  return from_matrix_cl<double>(x.matrix());
}

}  // namespace opencl
}  // namespace math
}  // namespace stan

#endif
#endif
