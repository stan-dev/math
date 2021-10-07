#ifndef STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_UNARY_HPP
#define STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_UNARY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_complex.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <utility>
#include <vector>

namespace stan {
namespace math {

/**
 * Base template class for vectorization of unary scalar functions
 * defined by a template class <code>F</code> to a scalar,
 * standard library vector, or Eigen dense matrix expression
 * template.
 *
 * <p>Specializations of this base define applications to scalars
 * (primitive or autodiff in the corresponding autodiff library
 * directories) or to standard library vectors of vectorizable
 * types (primitives, Eigen dense matrix expressions, or further
 * standard vectors).
 *
 * <p>Each specialization must define the typedef
 * <code>return_t</code> for the vectorized return type and the
 * function <code>apply</code> which defines the vectorization or
 * base application of the function defined statically by the
 * class F.  The function definition class F defines a static
 * function <code>fun()</code>, which defines the function's
 * behavior on scalars.
 *
 * @tparam F Type of function to apply.
 * @tparam T Type of argument to which function is applied.
 * @tparam ApplyZero If true, the function applied is assumed to return zero for
 *  inputs of zero and so sparse matrices will return sparse matrices. A value
 * of false will return a dense matrix for sparse matrices.
 */
template <typename F, typename T, bool ApplyZero = false,
          typename Enable = void>
struct apply_scalar_unary;

/**
 *
 * Template specialization for vectorized functions applying to
 * Eigen matrix arguments.
 *
 * @tparam F Type of function to apply.
 * @tparam T Type of argument to which function is applied.
 */
template <typename F, typename T, bool ApplyZero>
struct apply_scalar_unary<F, T, ApplyZero, require_eigen_t<T>> {
  /**
   * Type of underlying scalar for the matrix type T.
   */
  using scalar_t = value_type_t<T>;

  /**
   * Return the result of applying the function defined by the
   * template parameter F to the specified matrix argument.
   *
   * @tparam DenseMat A type derived from `Eigen::DenseBase`.
   * @param x Matrix to which operation is applied.
   * @return Componentwise application of the function specified
   * by F to the specified matrix.
   */
  template <typename DenseMat, require_eigen_dense_base_t<DenseMat>* = nullptr>
  static inline auto apply(const DenseMat& x) {
    return x.unaryExpr(
        [](scalar_t x) { return apply_scalar_unary<F, scalar_t>::apply(x); });
  }

  /**
   * Special case for `ApplyZero` set to true, returning a dense matrix. Return
   * the result of applying the function defined by the template parameter F to
   * the specified matrix argument.
   *
   * @param SparseMat A type derived from `Eigen::SparseMatrixBase`
   * @tparam NonZeroZero Shortcut trick for using class template for deduction,
   * should not be set manually.
   * @param x Matrix to which operation is applied.
   * @return Componentwise application of the function specified
   * by F to the specified matrix.
   */
  template <typename SparseMat, bool NonZeroZero = ApplyZero,
            require_t<bool_constant<NonZeroZero>>* = nullptr,
            require_eigen_sparse_base_t<SparseMat>* = nullptr>
  static inline auto apply(const SparseMat& x) {
    using val_t = value_type_t<SparseMat>;
    using triplet_t = Eigen::Triplet<val_t>;
    auto zeroed_val = apply_scalar_unary<F, scalar_t>::apply(val_t(0.0));
    const auto x_size = x.size();
    std::vector<triplet_t> triplet_list(x_size, triplet_t(0, 0, zeroed_val));
    for (Eigen::Index i = 0; i < x.rows(); ++i) {
      for (Eigen::Index j = 0; j < x.cols(); ++j) {
        // Column major order
        triplet_list[i * x.cols() + j] = triplet_t(i, j, zeroed_val);
      }
    }
    for (Eigen::Index k = 0; k < x.outerSize(); ++k) {
      for (typename SparseMat::InnerIterator it(x, k); it; ++it) {
        triplet_list[it.row() * x.cols() + it.col()] = triplet_t(
            it.row(), it.col(), apply_scalar_unary<F, scalar_t>::apply(it.value()));
      }
    }
    plain_type_t<SparseMat> ret(x.rows(), x.cols());
    ret.setFromTriplets(triplet_list.begin(), triplet_list.end());
    return ret;
  }

  /**
   * Special case for `ApplyZero` set to false, returning a sparse matrix.
   * Return the result of applying the function defined by the template
   * parameter F to the specified matrix argument.
   *
   * @tparam SparseMat A type derived from `Eigen::SparseMatrixBase`
   * @tparam NonZeroZero Shortcut trick for using class template for deduction,
   * should not be set manually.
   * @param x Matrix to which operation is applied.
   * @return Componentwise application of the function specified
   * by F to the specified matrix.
   */
  template <typename SparseMat, bool ReturnZeros = ApplyZero,
            require_t<bool_constant<!ReturnZeros>>* = nullptr,
            require_eigen_sparse_base_t<SparseMat>* = nullptr>
  static inline auto apply(const SparseMat& x) {
    auto ret = x.eval();
    for (Eigen::Index k = 0; k < x.outerSize(); ++k) {
      for (typename SparseMat::InnerIterator it(x, k), ret_it(ret, k); it;
           ++it, ++ret_it) {
        ret_it.valueRef() = apply_scalar_unary<F, scalar_t>::apply(it.value());
      }
    }
    return ret;
  }

  /**
   * Return type for applying the function elementwise to a matrix
   * expression template of type T.
   */
  using return_t = std::decay_t<decltype(
      apply_scalar_unary<F, T>::apply(std::declval<T>()))>;
};

/**
 * Template specialization for vectorized functions applying to
 * double arguments.
 *
 * @tparam F Type of function defining static apply function.
 */
template <typename F, typename T, bool ApplyZero>
struct apply_scalar_unary<F, T, ApplyZero, require_floating_point_t<T>> {
  /**
   * The return type, double.
   */
  using return_t = double;

  /**
   * Apply the function specified by F to the specified argument.
   * This is defined through a direct application of
   * <code>F::fun()</code>, which must be defined for double
   * arguments.
   *
   * @param x Argument scalar.
   * @return Result of applying F to the scalar.
   */
  static inline return_t apply(T x) { return F::fun(x); }
};

/**
 * Template specialization for vectorized functions applying to
 * complex arguments.
 *
 * @tparam F Type of function defining static apply function.
 */
template <typename F, typename T, bool ApplyZero>
struct apply_scalar_unary<F, T, ApplyZero, require_complex_t<T>> {
  /**
   * The return type, double.
   */
  using return_t = std::decay_t<T>;

  /**
   * Apply the function specified by F to the specified argument.
   * This is defined through a direct application of
   * <code>F::fun()</code>, which must be defined for double
   * arguments.
   *
   * @param x Argument scalar.
   * @return Result of applying F to the scalar.
   */
  static inline return_t apply(const T& x) { return F::fun(x); }
};

/**
 * Template specialization for vectorized functions applying to
 * integer arguments.  Although the argument is integer, the
 * return type is specified as double.  This allows promotion of
 * integers to doubles in vectorized functions, or in containers.
 *
 * @tparam F Type of function defining static apply function.
 */
template <typename F, typename T, bool ApplyZero>
struct apply_scalar_unary<F, T, ApplyZero, require_integral_t<T>> {
  /**
   * The return type, double.
   */
  using return_t = double;

  /**
   * Apply the function specified by F to the specified argument.
   * This is defined through a direct application of
   * <code>F::fun()</code>, which must be defined for double
   * arguments.
   *
   * @param x Argument scalar.
   * @return Result of applying F to the scalar.
   */
  static inline return_t apply(T x) { return F::fun(static_cast<double>(x)); }
};

/**
 * Template specialization for vectorized functions applying to
 * standard vector containers.  The lowest-level scalar type of
 * the argument will determine the return type.  Integers are
 * promoted to double values.
 *
 * @tparam F Class defining a static apply function.
 * @tparam T Type of element contained in standard vector.
 */
template <typename F, typename T, bool ApplyZero>
struct apply_scalar_unary<F, std::vector<T>, ApplyZero, void> {
  /**
   * Return type, which is calculated recursively as a standard
   * vector of the return type of the contained type T.
   */
  using return_t = typename std::vector<
      plain_type_t<typename apply_scalar_unary<F, T>::return_t>>;

  /**
   * Apply the function specified by F elementwise to the
   * specified argument.  This is defined recursively through this
   * class applied to elements of type T.
   *
   * @param x Argument container.
   * @return Elementwise application of F to the elements of the
   * container.
   */
  static inline return_t apply(const std::vector<T>& x) {
    return_t fx(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      fx[i] = apply_scalar_unary<F, T>::apply(x[i]);
    }
    return fx;
  }
};

}  // namespace math
}  // namespace stan
#endif
