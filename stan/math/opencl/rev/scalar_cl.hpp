#ifndef STAN_MATH_OPENCL_REV_SCALAR_CL_HPP
#define STAN_MATH_OPENCL_REV_SCALAR_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/callback_vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/scalar_cl.hpp>
#include <stan/math/opencl/scalar_cl_functions.hpp>
#include <stan/math/opencl/rev/arena_matrix_cl.hpp>
#include <stan/math/opencl/rev/vari.hpp>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
namespace opencl {

/** \ingroup opencl
 * A read-only device scalar whose buffer is owned by the autodiff arena, so it
 * can be captured by reverse-pass callbacks, whose destructors are not run.
 * The buffer is freed by `recover_memory()`. This is `arena_t` of the
 * primitive device scalars; constructing it copies the value on the device.
 */
template <>
class ScalarCl<arena_matrix_cl<double>>
    : public internal::scalar_cl_access<ScalarCl<arena_matrix_cl<double>>> {
 private:
  arena_matrix_cl<double> buf_;

 public:
  using value_type = double;

  /**
   * Copies the value of a primitive device scalar into the arena on the
   * device.
   * @tparam T type of the device scalar
   * @param x device scalar
   */
  template <typename T, require_prim_scalar_cl_t<T>* = nullptr,
            require_not_same_t<std::decay_t<T>,
                               ScalarCl<arena_matrix_cl<double>>>* = nullptr>
  ScalarCl(const T& x)  // NOLINT(runtime/explicit)
      : buf_(x.matrix()) {}

  ScalarCl(const ScalarCl& other) = default;
  ScalarCl& operator=(const ScalarCl& other) = default;

  /**
   * @return the backing 1x1 `matrix_cl`
   */
  inline const matrix_cl<double>& matrix() const noexcept {
    return static_cast<const matrix_cl<double>&>(buf_);
  }
};

/** \ingroup opencl
 * A var whose value and adjoint live on the OpenCL device.
 *
 * The autodiff node is a 1x1 `var_value<matrix_cl<double>>`, so the node and
 * its reverse-pass callbacks stay on the CPU while the value and adjoint stay
 * on the device. Copies share the autodiff node. Assignment and compound
 * assignment rebind to a new node; they never overwrite a value that earlier
 * operations may still use in the reverse pass.
 *
 * There are no implicit conversions to or from host values. Use the explicit
 * constructors to upload a value or a CPU `var`, and `opencl::to_host()` to
 * get a CPU `var`.
 */
template <>
class ScalarCl<var> {
 private:
  var_value<matrix_cl<double>> vi_;

 public:
  using value_type = var;

  /**
   * Constructs a constant device var holding zero.
   */
  ScalarCl() : vi_(std::move(ScalarCl<double>().matrix())) {}

  /**
   * Uploads a host value as a constant device var.
   * @param value value to upload
   */
  explicit ScalarCl(double value)
      : vi_(std::move(ScalarCl<double>(value).matrix())) {}

  /**
   * Constructs a constant device var from a device scalar, copying its value
   * on the device.
   * @tparam T type of the device scalar
   * @param value device scalar
   */
  template <typename T, require_prim_scalar_cl_t<T>* = nullptr>
  explicit ScalarCl(const T& value)
      : vi_(std::move(ScalarCl<double>(value).matrix())) {}

  /**
   * Moves a CPU var to the device.
   *
   * The value is uploaded now. In the reverse pass the device adjoint is read
   * back and added to the adjoint of `x`; this is the only transfer, and it is
   * made explicit by this constructor.
   * @param x CPU var
   */
  explicit ScalarCl(const var& x)
      : vi_(make_callback_var(std::move(ScalarCl<double>(x.val()).matrix()),
                              [x](vari_value<matrix_cl<double>>& res) mutable {
                                x.adj() += from_matrix_cl<double>(res.adj());
                              })) {}

  /**
   * Wraps an existing 1x1 autodiff node.
   * @param node autodiff node
   * @throw std::invalid_argument if the node is not 1x1
   */
  explicit ScalarCl(const var_value<matrix_cl<double>>& node) : vi_(node) {
    internal::check_scalar_cl_matrix(vi_.val());
  }

  ScalarCl(const ScalarCl& other) = default;
  ScalarCl(ScalarCl&& other) = default;
  ScalarCl& operator=(const ScalarCl& other) = default;
  ScalarCl& operator=(ScalarCl&& other) = default;

  /**
   * @return read-only view of the value, sharing its events
   */
  inline ScalarCl<const double&> val() const {
    return ScalarCl<const double&>(vi_.val());
  }

  /**
   * @return writable view of the adjoint, sharing its events
   */
  inline ScalarCl<double&> adj() const { return ScalarCl<double&>(vi_.adj()); }

  /**
   * @return the 1x1 autodiff node
   */
  inline const var_value<matrix_cl<double>>& node() const noexcept {
    return vi_;
  }

  /**
   * Adds to this device var by rebinding it to the sum.
   * @tparam T type of the argument
   * @param b value to add
   * @return this device var
   */
  template <typename T>
  inline ScalarCl& operator+=(T&& b);

  /**
   * Subtracts from this device var by rebinding it to the difference.
   * @tparam T type of the argument
   * @param b value to subtract
   * @return this device var
   */
  template <typename T>
  inline ScalarCl& operator-=(T&& b);

  /**
   * Multiplies this device var by rebinding it to the product.
   * @tparam T type of the argument
   * @param b value to multiply by
   * @return this device var
   */
  template <typename T>
  inline ScalarCl& operator*=(T&& b);

  /**
   * Divides this device var by rebinding it to the quotient.
   * @tparam T type of the argument
   * @param b value to divide by
   * @return this device var
   */
  template <typename T>
  inline ScalarCl& operator/=(T&& b);
};

/** \ingroup opencl
 * Returns a read-only view of the value of a device var. No data is
 * transferred.
 * @param x device var
 * @return view of the value
 */
inline ScalarCl<const double&> value_of(const ScalarCl<var>& x) {
  return x.val();
}

/** \ingroup opencl
 * Returns a writable view of the adjoint of a device var. No data is
 * transferred.
 * @param x device var
 * @return view of the adjoint
 */
inline ScalarCl<double&> adjoint_of(const ScalarCl<var>& x) { return x.adj(); }

/** \ingroup opencl
 * Moves a device var to a CPU var.
 *
 * The value is read back now, blocking until it is available. In the reverse
 * pass the CPU adjoint is added to the device adjoint without reading anything
 * back.
 * @param x device var
 * @return CPU var
 */
inline var to_host(const ScalarCl<var>& x) {
  return make_callback_var(
      to_host(x.val()), [x](const vari& res) mutable { x.adj() += res.adj(); });
}

/** \ingroup opencl
 * Creates a device var whose value is computed now and whose reverse pass is
 * given by a functor. The functor is called in the reverse pass with read-only
 * views of the adjoint and of the value of the result, and must add the
 * contributions of the result adjoint to the adjoints of the operands.
 *
 * The functor is stored in memory that is freed by `recover_memory()`, and
 * its destructor is called then, so it may capture device scalars and
 * `arena_matrix_cl` by value.
 * @tparam F type of the functor
 * @param value value of the result
 * @param functor reverse pass functor taking `(adjoint, value)`
 * @return device var
 */
template <typename F>
inline ScalarCl<var> make_callback_scalar_cl(ScalarCl<double>&& value,
                                             F&& functor) {
  return ScalarCl<var>(make_callback_var(
      std::move(value.matrix()),
      [functor
       = std::forward<F>(functor)](vari_value<matrix_cl<double>>& res) mutable {
        functor(ScalarCl<const double&>(res.adj()),
                ScalarCl<const double&>(res.val()));
      }));
}

namespace internal {

/**
 * Checks if a type can be an operand of an operation on device vars: a host
 * arithmetic value, a CPU var or a device scalar.
 */
template <typename T>
struct is_scalar_var
    : math::conjunction<is_var<std::decay_t<T>>,
                        std::is_floating_point<value_type_t<std::decay_t<T>>>> {
};

template <typename T>
struct is_scalar_cl_rev_operand
    : math::disjunction<std::is_arithmetic<std::decay_t<T>>, is_scalar_var<T>,
                        is_scalar_cl<T>> {};

/**
 * Checks if a type is a CPU var or a device var.
 */
template <typename T>
struct is_autodiff_scalar_cl_operand
    : math::disjunction<is_rev_scalar_cl<T>, is_scalar_var<T>> {};

/**
 * Enables a template if all types are host arithmetic values, CPU vars or
 * device scalars, at least one is a device scalar and at least one is
 * autodiff.
 */
template <typename... Types>
using require_scalar_cl_rev_operands_t = require_t<math::conjunction<
    is_scalar_cl_rev_operand<Types>...,
    math::disjunction<is_scalar_cl<Types>...>,
    math::disjunction<is_autodiff_scalar_cl_operand<Types>...>>>;

/**
 * Prepares an operand for use in the forward and reverse pass of an operation
 * on device vars. CPU vars are moved to the device, device scalars are copied
 * into the arena so the reverse pass sees the value used in the forward pass,
 * and device vars and host values are kept as they are.
 */
inline ScalarCl<var> to_rev_operand(const ScalarCl<var>& x) { return x; }
inline ScalarCl<var> to_rev_operand(const var& x) { return ScalarCl<var>(x); }
inline double to_rev_operand(double x) { return x; }
template <typename T, require_prim_scalar_cl_t<T>* = nullptr>
inline arena_matrix_cl<double> to_rev_operand(const T& x) {
  return arena_matrix_cl<double>(x.matrix());
}

/**
 * Returns the value of a prepared operand as something the kernel generator
 * can use: a view of a device value or a host value.
 */
inline ScalarCl<const double&> operand_value(const ScalarCl<var>& x) {
  return x.val();
}
inline ScalarCl<const double&> operand_value(const arena_matrix_cl<double>& x) {
  return ScalarCl<const double&>(static_cast<const matrix_cl<double>&>(x));
}
inline double operand_value(double x) { return x; }

/**
 * Adds a derivative expression to the adjoint of a prepared operand if it is
 * a device var. Other operands have no adjoint.
 */
template <typename Expr>
inline void add_to_adjoint(const ScalarCl<var>& x, Expr&& expr) {
  x.adj() += std::forward<Expr>(expr);
}
template <typename T, typename Expr, require_not_rev_scalar_cl_t<T>* = nullptr>
inline void add_to_adjoint(const T&, Expr&&) {}

/**
 * Builds a device var from a binary operation.
 * @param a first operand
 * @param b second operand
 * @param value_f computes the value from the values of the operands
 * @param da_f computes the derivative expression for `a` from the result
 * adjoint, the result value and the operand values
 * @param db_f computes the derivative expression for `b`
 * @return device var
 */
template <typename T_a, typename T_b, typename F_val, typename F_da,
          typename F_db>
inline ScalarCl<var> scalar_cl_binary(const T_a& a, const T_b& b,
                                      F_val&& value_f, F_da&& da_f,
                                      F_db&& db_f) {
  auto a_op = to_rev_operand(a);
  auto b_op = to_rev_operand(b);
  return make_callback_scalar_cl(
      value_f(operand_value(a_op), operand_value(b_op)),
      [a_op, b_op, da_f, db_f](const auto& res_adj,
                               const auto& res_val) mutable {
        auto a_val = operand_value(a_op);
        auto b_val = operand_value(b_op);
        add_to_adjoint(a_op, da_f(res_adj, res_val, a_val, b_val));
        add_to_adjoint(b_op, db_f(res_adj, res_val, a_val, b_val));
      });
}

/**
 * Builds a device var from a unary operation.
 * @param a device var operand
 * @param value_f computes the value from the value of the operand
 * @param da_f computes the derivative expression from the result adjoint, the
 * result value and the operand value
 * @return device var
 */
template <typename F_val, typename F_da>
inline ScalarCl<var> scalar_cl_unary(const ScalarCl<var>& a, F_val&& value_f,
                                     F_da&& da_f) {
  return make_callback_scalar_cl(
      value_f(a.val()),
      [a, da_f](const auto& res_adj, const auto& res_val) mutable {
        a.adj() += da_f(res_adj, res_val, a.val());
      });
}

}  // namespace internal

/** \ingroup opencl
 * Adds scalars, at least one of which is a device scalar and at least one of
 * which is autodiff.
 * @return device var holding the sum
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_rev_operands_t<T_a, T_b>* = nullptr>
inline ScalarCl<var> operator+(const T_a& a, const T_b& b) {
  return internal::scalar_cl_binary(
      a, b, [](const auto& x, const auto& y) { return x + y; },
      [](const auto& g, const auto&, const auto&, const auto&) {
        return as_operation_cl(g);
      },
      [](const auto& g, const auto&, const auto&, const auto&) {
        return as_operation_cl(g);
      });
}

/** \ingroup opencl
 * Subtracts scalars, at least one of which is a device scalar and at least one
 * of which is autodiff.
 * @return device var holding the difference
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_rev_operands_t<T_a, T_b>* = nullptr>
inline ScalarCl<var> operator-(const T_a& a, const T_b& b) {
  return internal::scalar_cl_binary(
      a, b, [](const auto& x, const auto& y) { return x - y; },
      [](const auto& g, const auto&, const auto&, const auto&) {
        return as_operation_cl(g);
      },
      [](const auto& g, const auto&, const auto&, const auto&) {
        return -as_operation_cl(g);
      });
}

/** \ingroup opencl
 * Multiplies scalars, at least one of which is a device scalar and at least one
 * of which is autodiff.
 * @return device var holding the product
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_rev_operands_t<T_a, T_b>* = nullptr>
inline ScalarCl<var> operator*(const T_a& a, const T_b& b) {
  return internal::scalar_cl_binary(
      a, b, [](const auto& x, const auto& y) { return x * y; },
      [](const auto& g, const auto&, const auto&, const auto& y) {
        return elt_multiply(as_operation_cl(g), as_operation_cl(y));
      },
      [](const auto& g, const auto&, const auto& x, const auto&) {
        return elt_multiply(as_operation_cl(g), as_operation_cl(x));
      });
}

/** \ingroup opencl
 * Divides scalars, at least one of which is a device scalar and at least one
 * of which is autodiff.
 * @return device var holding the quotient
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_rev_operands_t<T_a, T_b>* = nullptr>
inline ScalarCl<var> operator/(const T_a& a, const T_b& b) {
  return internal::scalar_cl_binary(
      a, b, [](const auto& x, const auto& y) { return x / y; },
      [](const auto& g, const auto&, const auto&, const auto& y) {
        return elt_divide(as_operation_cl(g), as_operation_cl(y));
      },
      [](const auto& g, const auto& r, const auto&, const auto& y) {
        return -elt_divide(elt_multiply(as_operation_cl(g), as_operation_cl(r)),
                           as_operation_cl(y));
      });
}

/** \ingroup opencl
 * Negates a device var.
 * @param a device var
 * @return device var holding the negation
 */
inline ScalarCl<var> operator-(const ScalarCl<var>& a) {
  return internal::scalar_cl_unary(
      a, [](const auto& x) { return -x; },
      [](const auto& g, const auto&, const auto&) {
        return -as_operation_cl(g);
      });
}

/** \ingroup opencl
 * Exponential of a device var.
 * @param a device var
 * @return device var holding the exponential
 */
inline ScalarCl<var> exp(const ScalarCl<var>& a) {
  return internal::scalar_cl_unary(
      a, [](const auto& x) { return exp(x); },
      [](const auto& g, const auto& r, const auto&) {
        return elt_multiply(as_operation_cl(g), as_operation_cl(r));
      });
}

/** \ingroup opencl
 * Natural logarithm of a device var.
 * @param a device var
 * @return device var holding the logarithm
 */
inline ScalarCl<var> log(const ScalarCl<var>& a) {
  return internal::scalar_cl_unary(
      a, [](const auto& x) { return log(x); },
      [](const auto& g, const auto&, const auto& x) {
        return elt_divide(as_operation_cl(g), as_operation_cl(x));
      });
}

/** \ingroup opencl
 * Square root of a device var.
 * @param a device var
 * @return device var holding the square root
 */
inline ScalarCl<var> sqrt(const ScalarCl<var>& a) {
  return internal::scalar_cl_unary(
      a, [](const auto& x) { return sqrt(x); },
      [](const auto& g, const auto& r, const auto&) {
        return elt_divide(0.5 * as_operation_cl(g), as_operation_cl(r));
      });
}

/** \ingroup opencl
 * Square of a device var.
 * @param a device var
 * @return device var holding the square
 */
inline ScalarCl<var> square(const ScalarCl<var>& a) {
  return internal::scalar_cl_unary(
      a, [](const auto& x) { return square(x); },
      [](const auto& g, const auto&, const auto& x) {
        return 2.0 * elt_multiply(as_operation_cl(g), as_operation_cl(x));
      });
}

namespace internal {
/**
 * Enables a template if one type is a device scalar, the other is a matrix
 * expression, and at least one of them is autodiff.
 */
template <typename T_a, typename T_b>
using require_scalar_cl_and_matrix_rev_t = require_t<math::conjunction<
    math::disjunction<
        math::conjunction<is_scalar_cl<T_a>,
                          is_nonscalar_prim_or_rev_kernel_expression<T_b>>,
        math::conjunction<is_nonscalar_prim_or_rev_kernel_expression<T_a>,
                          is_scalar_cl<T_b>>>,
    math::disjunction<is_var<scalar_type_t<T_a>>, is_var<scalar_type_t<T_b>>>>>;
}  // namespace internal

/** \ingroup opencl
 * Adds a device scalar and a matrix, at least one of which is autodiff.
 * @return var matrix holding the sum
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_and_matrix_rev_t<T_a, T_b>* = nullptr>
inline auto operator+(const T_a& a, const T_b& b) {
  return add(a, b);
}

/** \ingroup opencl
 * Subtracts a device scalar and a matrix, at least one of which is autodiff.
 * @return var matrix holding the difference
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_and_matrix_rev_t<T_a, T_b>* = nullptr>
inline auto operator-(const T_a& a, const T_b& b) {
  return subtract(a, b);
}

/** \ingroup opencl
 * Multiplies a device scalar and a matrix, at least one of which is autodiff.
 * @return var matrix holding the product
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_and_matrix_rev_t<T_a, T_b>* = nullptr>
inline auto operator*(const T_a& a, const T_b& b) {
  return multiply(a, b);
}

/** \ingroup opencl
 * Divides a device scalar and a matrix element-wise, at least one of which is
 * autodiff.
 * @return var matrix holding the quotient
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_and_matrix_rev_t<T_a, T_b>* = nullptr>
inline auto operator/(const T_a& a, const T_b& b) {
  if constexpr (is_scalar_cl<T_b>::value) {
    return divide(a, b);
  } else {
    return elt_divide(a, b);
  }
}

template <typename T>
inline ScalarCl<var>& ScalarCl<var>::operator+=(T&& b) {
  *this = *this + b;
  return *this;
}

template <typename T>
inline ScalarCl<var>& ScalarCl<var>::operator-=(T&& b) {
  *this = *this - b;
  return *this;
}

template <typename T>
inline ScalarCl<var>& ScalarCl<var>::operator*=(T&& b) {
  *this = *this * b;
  return *this;
}

template <typename T>
inline ScalarCl<var>& ScalarCl<var>::operator/=(T&& b) {
  *this = *this / b;
  return *this;
}

}  // namespace opencl
}  // namespace math
}  // namespace stan

#endif
#endif
