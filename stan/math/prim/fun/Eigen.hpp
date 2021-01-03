#ifndef STAN_MATH_PRIM_FUN_EIGEN_HPP
#define STAN_MATH_PRIM_FUN_EIGEN_HPP

namespace Eigen {

namespace internal {
// Workaround for SFINAE regression with GCC 8.3
template <typename>
struct has_member {
  typedef int type;
};

/**
 * Reimplements is_fvar without requiring external math headers
 *
 * decltype((void)(T::d_)) is a pre C++17 replacement for
 * std::void_t<decltype(T::d_)>
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template <typename T,
          typename has_member<decltype(std::decay_t<T>::d_)>::type = 0>
struct is_fvar : std::true_type {};

/**
 * Reimplements is_var without requiring external math headers
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template <typename T,
          typename has_member<decltype(std::decay_t<T>::vi_)>::type = 0>
struct is_var : std::true_type {};

/**
 * Reimplements is_vari without requiring external math headers
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template <typename T, typename has_member<decltype(
                          std::remove_pointer_t<std::decay_t<T>>::adj_)>::type
                      = 0>
struct is_vari : std::true_type {};

struct val_Op;
struct val_ref_Op;

struct d_Op;
struct d_ref_Op;

struct adj_Op;
struct adj_ref_Op;

struct vi_Op;
struct vi_ref_Op;

}  // namespace internal
}  // namespace Eigen
#ifdef EIGEN_MATRIXBASE_PLUGIN
#ifndef EIGEN_STAN_MATRIXBASE_PLUGIN
#error "Stan uses Eigen's EIGEN_MATRIXBASE_PLUGIN macro. To use your own "
"plugin add the eigen_plugin.h file to your plugin."
#endif
#else
#define EIGEN_MATRIXBASE_PLUGIN "stan/math/prim/eigen_plugins.h"
#endif

#ifdef EIGEN_ARRAYBASE_PLUGIN
#ifndef EIGEN_STAN_ARRAYBASE_PLUGIN
#error "Stan uses Eigen's EIGEN_ARRAYBASE_PLUGIN macro. To use your own "
    "plugin add the eigen_plugin.h file to your plugin."
#endif
#else
#define EIGEN_ARRAYBASE_PLUGIN "stan/math/prim/eigen_plugins.h"
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/src/Core/NumTraits.h>

    namespace Eigen {

  namespace internal {

  /**
   * Structure to return a view to the values in a var, vari*, and fvar<T>.
   * To identify the correct member to call for a given input, templates
   * check a combination of whether the input is a pointer (i.e. vari*)
   * and/or whether the input has member ".d_" (i.e. fvar).
   *
   * There are two methods for returning doubles unchanged. One which takes a
   * reference to a double and returns the same reference, used when 'chaining'
   * methods (i.e. A.adj().val()). The other for passing and returning by value,
   * used directly with matrices of doubles (i.e. A.val(), where A is of type
   * MatrixXd).
   *
   * For definitions of EIGEN_EMPTY_STRUCT_CTOR, EIGEN_DEVICE_FUNC, and
   * EIGEN_STRONG_INLINE; see:
   * https://eigen.tuxfamily.org/dox/XprHelper_8h_source.html
   */
  struct val_Op {
    EIGEN_EMPTY_STRUCT_CTOR(val_Op);

    // Returns value from a vari*
    template <typename T, std::enable_if_t<is_vari<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(const T& v) const {
      return numext::real(v->val_);
    }

    // Returns value from a var
    template <typename T, std::enable_if_t<is_var<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(const T& v) const {
      return numext::real(v.vi_->val_);
    }

    // Returns value from an fvar
    template <typename T, std::enable_if_t<is_fvar<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(const T& v) const {
      return v.val_;
    }

    // Returns double unchanged from input (by reference)
    template <
        typename T,
        std::enable_if_t<std::is_arithmetic<std::decay_t<T>>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(const T& v) const {
      return numext::real(v);
    }
  };
  template <>
  struct functor_traits<val_Op> {
    enum { Cost = 0, PacketAccess = false };
  };

  struct val_ref_Op {
    EIGEN_EMPTY_STRUCT_CTOR(val_ref_Op);

    // Returns value from a vari*
    template <typename T, std::enable_if_t<is_vari<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto& operator()(const T& v) const {
      return numext::real_ref(*const_cast<double*>(&(v->val_)));
    }

    // Returns value from a var
    template <typename T, std::enable_if_t<is_var<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto& operator()(const T& v) const {
      return numext::real_ref(*const_cast<double*>(&(v.vi_->val_)));
    }

    // Returns value from an fvar
    template <typename T, std::enable_if_t<is_fvar<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto& operator()(const T& v) const {
      return const_cast<typename std::decay_t<T>::Scalar&>(v.val_);
    }

    // Returns double unchanged from input (by reference)
    template <
        typename T,
        std::enable_if_t<std::is_arithmetic<std::decay_t<T>>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto& operator()(const T& v) const {
      return numext::real_ref(*const_cast<double*>(&v));
    }
  };
  template <>
  struct functor_traits<val_ref_Op> {
    enum { Cost = 0, PacketAccess = false };
  };

  /**
   * Structure to return tangent from an fvar.
   */
  struct d_Op {
    EIGEN_EMPTY_STRUCT_CTOR(d_Op);

    // Returns tangent from an fvar
    template <typename T, std::enable_if_t<is_fvar<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(const T& v) const {
      return v.d_;
    }
  };
  template <>
  struct functor_traits<d_Op> {
    enum { Cost = 0, PacketAccess = false };
  };

  struct d_ref_Op {
    EIGEN_EMPTY_STRUCT_CTOR(d_ref_Op);

    // Returns tangent from an fvar
    template <typename T, std::enable_if_t<is_fvar<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto& operator()(const T& v) const {
      return const_cast<typename std::decay_t<T>::Scalar&>(v.d_);
    }
  };
  template <>
  struct functor_traits<d_ref_Op> {
    enum { Cost = 0, PacketAccess = false };
  };

  /**
   * Structure to return adjoints from var and vari*.
   */
  struct adj_Op {
    EIGEN_EMPTY_STRUCT_CTOR(adj_Op);

    // Returns adjoint from a vari*
    template <typename T, std::enable_if_t<is_vari<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(const T& v) const {
      return v->adj_;
    }

    // Returns adjoint from a var
    template <typename T, std::enable_if_t<is_var<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(const T& v) const {
      return v.vi_->adj_;
    }
  };
  template <>
  struct functor_traits<adj_Op> {
    enum { Cost = 0, PacketAccess = false };
  };

  struct adj_ref_Op {
    EIGEN_EMPTY_STRUCT_CTOR(adj_ref_Op);

    // Returns adjoint from a vari*
    template <typename T, std::enable_if_t<is_vari<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto& operator()(const T& v) const {
      return numext::real_ref(*const_cast<double*>(&(v->adj_)));
    }

    // Returns adjoint from a var
    template <typename T, std::enable_if_t<is_var<T>::value>* = nullptr>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto& operator()(const T& v) const {
      return numext::real_ref(*const_cast<double*>(&(v.vi_->adj_)));
    }
  };
  template <>
  struct functor_traits<adj_ref_Op> {
    enum { Cost = 0, PacketAccess = false };
  };

  /**
   * Structure to return vari* from a var.
   */
  struct vi_Op {
    EIGEN_EMPTY_STRUCT_CTOR(vi_Op);

    // Returns vari* from a var
    template <typename T>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(const T& v) const {
      return v.vi_;
    }
  };
  template <>
  struct functor_traits<vi_Op> {
    enum { Cost = 0, PacketAccess = false };
  };

  /**
   * Structure to return vari* from a var.
   */
  struct vi_ref_Op {
    EIGEN_EMPTY_STRUCT_CTOR(vi_ref_Op);

    // Returns vari* from a var
    template <typename T>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto& operator()(const T& v) const {
      return const_cast<std::decay_t<decltype(v.vi_)>&>(v.vi_);
    }
  };
  template <>
  struct functor_traits<vi_ref_Op> {
    enum { Cost = 0, PacketAccess = false };
  };

  }  // namespace internal

  /**
   * Traits specialization for Eigen binary operations for `int`
   * and `double` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<int, double, BinaryOp> {
    using ReturnType = double;
  };

  /**
   * Traits specialization for Eigen binary operations for `double`
   * and `int` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<double, int, BinaryOp> {
    using ReturnType = double;
  };

  /**
   * Traits specialization for Eigen binary operations for `int`
   * and complex `double` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<int, std::complex<double>, BinaryOp> {
    using ReturnType = std::complex<double>;
  };

  /**
   * Traits specialization for Eigen binary operations for complex
   * `double` and `int` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<std::complex<double>, int, BinaryOp> {
    using ReturnType = std::complex<double>;
  };

}  // namespace Eigen


#endif
