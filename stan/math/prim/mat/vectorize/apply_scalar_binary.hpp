#ifndef STAN_MATH_PRIM_MAT_VECTORIZE_APPLY_BINARY_SCALAR_HPP
#define STAN_MATH_PRIM_MAT_VECTORIZE_APPLY_BINARY_SCALAR_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/err/check_matching_dims.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <Eigen/Core>
#include <vector>

namespace stan {

namespace math {

/**
 * Base template class for vectorization of binary scalar functions
 * defined by a template class <code>F</code> to scalars,
 * standard library vectors, or Eigen dense matrix expression
 * templates.
 *
 * <p>The base class applies to any to scalars
 * (primitive or autodiff in the corresponding autodiff library
 * directories). Specializations define applications to dense matrix
 * expression template or to standard library vectors of vectorizable
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
 * @tparam T1 Type of first argument to which function is applied.
 * @tparam T2 Type of second argument to which function is applied.
 */
template <typename F, typename T1, typename T2>
struct apply_scalar_binary {
  /**
   * Return type for applying the function to the scalars
   * of types T1 and T2.
   */
  typedef typename return_type<T1, T2>::type return_t;

  /**
   * Return the result of applying the function defined by the
   * template parameter F to the specified scalars.
   *
   * @param x Scalar to which operation is applied.
   * @param y Scalar to which operation is applied.
   * @return Application of the function specified
   * by F to the specified scalars.
   */
  static inline return_t apply(const T1& x, const T2& y) {
    return F::fun(x, y);
  }
};

template <typename F, typename T1, typename T2, int R1, int C1, int R2, int C2>
struct apply_scalar_binary<F, Eigen::Matrix<T1, R1, C1>,
                           Eigen::Matrix<T2, R2, C2> > {
  /**
   * Return type, which is calculated recursively as a matrix
   * of the return type of the contained types T1 and T2.
   */
  typedef Eigen::Matrix<typename apply_scalar_binary<F, T1, T2>::return_t, R1,
                        C1>
      return_t;

  /**
   * Return the result of applying the function defined by the
   * template parameter F to the specified matrix arguments.
   *
   * @param x Matrix to which operation is applied.
   * @param y Matrix to which operation is applied.
   * @return Componentwise application of the function specified
   * by F to the specified matrices.
   */
  static inline return_t apply(const Eigen::Matrix<T1, R1, C1>& x,
                               const Eigen::Matrix<T2, R2, C2>& y) {
    check_matching_dims<true>("binary vectorization", "x", x, "y", y);
    return_t result(x.rows(), x.cols());
    for (int i = 0; i < x.size(); ++i)
      result(i) = apply_scalar_binary<F, T1, T2>::apply(x(i), y(i));
    return result;
  }
};

template <typename F, typename T1, typename T2, int R1, int C1, int R2, int C2>
struct apply_scalar_binary<F, Eigen::Block<Eigen::Matrix<T1, R1, C1> >,
                           Eigen::Block<Eigen::Matrix<T2, R2, C2> > > {
  /**
   * Return type, which is calculated recursively as an Eigen::Matrix
   * of the return type of the contained types T1 and T2.
   */
  typedef Eigen::Matrix<typename apply_scalar_binary<F, T1, T2>::return_t, R1,
                        C1>
      return_t;

  /**
   * Apply the function specified by F elementwise to the
   * specified argument. This is defined recursively through this
   * class by passing the arguments as matrices to the
   * <code>Eigen::Matrix</code> specialization.
   *
   * @param x <code>Eigen::Block</code> argument container.
   * @param y <code>Eigen::Block argument</code> container.
   * @return Elementwise application of F to the elements of the
   * containers.
   */
  static inline return_t apply(
      const Eigen::Block<Eigen::Matrix<T1, R1, C1> >& x,
      const Eigen::Block<Eigen::Matrix<T2, R2, C2> >& y) {
    check_matching_dims<true>(
        "binary vectorization", "x",
        static_cast<Eigen::Matrix<T1, R1, C1> >(x.eval()), "y",
        static_cast<Eigen::Matrix<T1, R1, C1> >(x.eval()));
    return_t result(x.rows(), x.cols());
    for (int i = 0; i < x.size(); ++i)
      result(i)
          = apply_scalar_binary<F, T1, T2>::apply(x.eval()(i), y.eval()(i));
    return result;
  }
};

template <typename F, typename T1, typename T2>
struct apply_scalar_binary<F, std::vector<T1>, std::vector<T2> > {
  /**
   * Return type, which is calculated recursively as a standard
   * vector of the return type of the contained types T1 and T2.
   */
  typedef
      typename std::vector<typename apply_scalar_binary<F, T1, T2>::return_t>
          return_t;

  /**
   * Apply the function specified by F elementwise to the
   * specified arguments.  This is defined recursively through this
   * class applied to elements of types T1 and T2.
   *
   * @param x Argument container.
   * @param y Argument container.
   * @return Elementwise application of F to the elements of the
   * containers.
   */
  static inline return_t apply(const std::vector<T1>& x,
                               const std::vector<T2>& y) {
    using stan::math::check_size_match;

    check_size_match("binary vectorization", "x", x.size(), "y", y.size());
    return_t fx(x.size());
    for (size_t i = 0; i < x.size(); ++i)
      fx[i] = apply_scalar_binary<F, T1, T2>::apply(x[i], y[i]);
    return fx;
  }
};
}  // namespace math
}  // namespace stan
#endif
