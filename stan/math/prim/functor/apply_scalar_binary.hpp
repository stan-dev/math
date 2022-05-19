#ifndef STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_BINARY_HPP
#define STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_BINARY_HPP

#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/num_elements.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Base template function for vectorization of binary scalar functions
 * defined by applying a functor to a combination of scalars,
 * containers of matching sizes, or a combination of a scalar and a container.
 * These containers can be a standard library vector, Eigen dense
 * matrix expression template, or container of these. For each specialisation,
 * the same type as the largest (dimension) input is returned.
 *
 * This base template function takes (and returns) scalars.
 *
 * @tparam Scalar1 Type of first argument to which functor is applied.
 * @tparam Scalar2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First input to which operation is applied.
 * @param y Second input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Scalar with result of applying functor to input.
 */
template <typename Scalar1, typename Scalar2, typename F,
          require_all_stan_scalar_t<Scalar1, Scalar2>* = nullptr>
inline auto apply_scalar_binary(const Scalar1& x, const Scalar2& y, const F& f) {
  return f(x, y);
}

/**
 * Specialisation for use with two Eigen inputs. Eigen's binaryExpr framework
 * is used for more efficient indexing of both row- and column-major inputs
 * without separate loops.
 *
 * @tparam Eig1 Type of first argument to which functor is applied.
 * @tparam Eig2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First Eigen input to which operation is applied.
 * @param y Second Eigen input to which operation is applied.
 * @param f functor to apply to Eigen input.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename Eig1, typename Eig2, typename F,
          require_all_eigen_t<Eig1, Eig2>* = nullptr>
inline auto apply_scalar_binary(const Eig1& x, const Eig2& y, const F& f) {
  check_matching_dims("Binary function", "x", x, "y", y);
  return x.binaryExpr(y, f);
}

/**
 * Specialisation for use with one Eigen vector (row or column) and
 * a one-dimensional std::vector of integer types
 *
 * @tparam EigVec Type of first argument to which functor is applied.
 * @tparam StdVec Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Eigen input to which operation is applied.
 * @param y Integer std::vector input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename EigVec, typename StdVec, typename F,
          require_eigen_vector_vt<is_stan_scalar, EigVec>* = nullptr,
          require_std_vector_vt<std::is_integral, StdVec>* = nullptr>
inline auto apply_scalar_binary(const EigVec& x, const StdVec& y, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  using int_vec_t = promote_scalar_t<value_type_t<StdVec>, plain_type_t<EigVec>>;
  Eigen::Map<const int_vec_t> y_map(y.data(), y.size());
  return x.binaryExpr(y_map, f);
}

/**
 * Specialisation for use with a one-dimensional std::vector of integer types
 * and one Eigen vector (row or column).
 *
 * @tparam StdVec Type of first argument to which functor is applied.
 * @tparam StdVec Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Integer std::vector input to which operation is applied.
 * @param y Eigen input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename StdVec, typename StdVec, typename F,
          require_std_vector_vt<std::is_integral, StdVec>* = nullptr,
          require_eigen_vector_vt<is_stan_scalar, StdVec>* = nullptr>
inline auto apply_scalar_binary(const StdVec& x, const StdVec& y, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  using int_vec_t = promote_scalar_t<value_type_t<StdVec>, plain_type_t<StdVec>>;
  Eigen::Map<const int_vec_t> x_map(x.data(), x.size());
  return x_map.binaryExpr(y, f);
}

/**
 * Specialisation for use with one Eigen matrix and
 * a two-dimensional std::vector of integer types
 *
 * @tparam EigMat Type of first argument to which functor is applied.
 * @tparam Container Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Eigen matrix input to which operation is applied.
 * @param y Nested integer std::vector input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename EigMat, typename Container, typename F,
          require_eigen_matrix_dynamic_vt<is_stan_scalar, EigMat>* = nullptr,
          require_std_vector_vt<is_std_vector, Container>* = nullptr,
          require_std_vector_st<std::is_integral, Container>* = nullptr>
inline auto apply_scalar_binary(const EigMat& x, const Container& y, const F& f) {
  if (num_elements(x) != num_elements(y)) {
    std::ostringstream msg;
    msg << "Inputs to vectorized binary function must match in"
        << " size if one is not a scalar";
    throw std::invalid_argument(msg.str());
  }
  using return_st = decltype(f(x(0), y[0][0]));
  Eigen::Matrix<return_st, Eigen::Dynamic, Eigen::Dynamic> result(x.rows(),
                                                                  x.cols());
  for (size_t i = 0; i < y.size(); ++i) {
    result.row(i) = apply_scalar_binary(x.row(i).transpose(),
                                        as_column_vector_or_scalar(y[i]), f);
  }
  return result;
}

/**
 * Specialisation for use with a two-dimensional std::vector of integer types
 * and one Eigen matrix.
 *
 * @tparam StdVecVec Type of first argument to which functor is applied.
 * @tparam EigMat Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Nested integer std::vector input to which operation is applied.
 * @param y Eigen matrix input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename StdVecVec, typename EigMat, typename F,
          require_std_vector_vt<is_std_vector, StdVecVec>* = nullptr,
          require_std_vector_st<std::is_integral, StdVecVec>* = nullptr,
          require_eigen_matrix_dynamic_vt<is_stan_scalar, EigMat>* = nullptr>
inline auto apply_scalar_binary(const StdVecVec& x, const EigMat& y, const F& f) {
  if (num_elements(x) != num_elements(y)) {
    std::ostringstream msg;
    msg << "Inputs to vectorized binary function must match in"
        << " size if one is not a scalar";
    throw std::invalid_argument(msg.str());
  }
  using return_st = decltype(f(x[0][0], y(0)));
  Eigen::Matrix<return_st, Eigen::Dynamic, Eigen::Dynamic> result(y.rows(),
                                                                  y.cols());
  for (size_t i = 0; i < x.size(); ++i) {
    result.row(i) = apply_scalar_binary(as_column_vector_or_scalar(x[i]),
                                        y.row(i).transpose(), f);
  }
  return result;
}

/**
 * Specialisation for use when the first input is an Eigen type and the second
 * is a scalar. Eigen's unaryExpr framework is used for more efficient indexing
 * of both row- and column-major inputs. The unaryExpr framework also allows
 * for the scalar to be captured and applied to each element in the Eigen
 * object.
 *
 * @tparam Eig Type of Eigen object to which functor is applied.
 * @tparam Scalar Type of scalar to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Eigen input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @param f functor to apply to Eigen and scalar inputs.
 * @return Eigen object with result of applying functor to inputs.
 *
 * Note: The return expresssion needs to be evaluated, otherwise the captured
 *         function and scalar fall out of scope.
 */
template <typename Eig, typename Scalar, typename F, require_eigen_t<Eig>* = nullptr,
          require_stan_scalar_t<Scalar>* = nullptr>
inline auto apply_scalar_binary(const Eig& x, const Scalar& y, const F& f) {
  return x.unaryExpr([&f, y](const auto& v) { return f(v, y); });
}

/**
 * Specialisation for use when the first input is an scalar and the second is
 * an Eigen type. Eigen's unaryExpr framework is used for more efficient
 * indexing of both row- and column-major inputs. The unaryExpr framework also
 * allows for the scalar to be captured and applied to each element in the
 * Eigen object.
 *
 * @tparam Scalar Type of scalar to which functor is applied.
 * @tparam Eig Type of Eigen object to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Scalar input to which operation is applied.
 * @param y Eigen input to which operation is applied.
 * @param f functor to apply to Eigen and scalar inputs.
 * @return Eigen object with result of applying functor to inputs.
 *
 * Note: The return expresssion needs to be evaluated, otherwise the captured
 *         function and scalar fall out of scope.
 */
template <typename Scalar, typename Eig, typename F,
          require_stan_scalar_t<Scalar>* = nullptr, require_eigen_t<Eig>* = nullptr>
inline auto apply_scalar_binary(const Scalar& x, const Eig& y, const F& f) {
  return y.unaryExpr([&f, x](const auto& v) { return f(x, v); });
}

/**
 * Specialisation for use with (non-nested) std::vectors. Inputs are mapped
 * to Eigen column vectors and then the result is evaluated directly into the
 * returned std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam StdVec1 Type of first std::vector to which functor is applied.
 * @tparam StdVec2 Type of second std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First std::vector input to which operation is applied.
 * @param y Second std::vector input to which operation is applied.
 * @param f functor to apply to std::vector inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename StdVec1, typename StdVec2, typename F,
          require_all_std_vector_vt<is_stan_scalar, StdVec1, StdVec2>* = nullptr>
inline auto apply_scalar_binary(const StdVec1& x, const StdVec2& y, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  decltype(auto) x_vec = as_column_vector_or_scalar(x);
  decltype(auto) y_vec = as_column_vector_or_scalar(y);
  using T_return = std::decay_t<decltype(f(x[0], y[0]))>;
  std::vector<T_return> result(x.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = x_vec.binaryExpr(y_vec, f);
  return result;
}

/**
 * Specialisation for use when the first input is a (non-nested) std::vector
 * and the second is a scalar. The std::vector input is mapped to an Eigen
 * column vector and then the result is evaluated directly into the returned
 * std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam StdVecScalar Type of std::vector to which functor is applied.
 * @tparam Scalar Type of scalar to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x std::vector input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @param f functor to apply to std::vector and scalar inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename StdVecScalar, typename Scalar, typename F,
          require_std_vector_vt<is_stan_scalar, StdVecScalar>* = nullptr,
          require_stan_scalar_t<Scalar>* = nullptr>
inline auto apply_scalar_binary(const StdVecScalar& x, const StdVec& y, const F& f) {
  decltype(auto) x_vec = as_column_vector_or_scalar(x);
  using T_return = std::decay_t<decltype(f(x[0], y))>;
  std::vector<T_return> result(x.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = x_vec.unaryExpr([&f, &y](const auto& v) { return f(v, y); });
  return result;
}

/**
 * Specialisation for use when the first input is a scalar and the second is a
 * (non-nested) std::vector. The std::vector input is mapped to an Eigen
 * column vector and then the result is evaluated directly into the returned
 * std::vector (via Eigen::Map).
 *
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam Scalar Type of scalar to which functor is applied.
 * @tparam StdVec Type of std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Scalar input to which operation is applied.
 * @param y std::vector input to which operation is applied.
 * @param f functor to apply to std::vector and scalar inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename Scalar, typename StdVec, typename F,
          require_stan_scalar_t<Scalar>* = nullptr,
          require_std_vector_vt<is_stan_scalar, StdVec>* = nullptr>
inline auto apply_scalar_binary(const Scalar& x, const StdVec& y, const F& f) {
  decltype(auto) y_vec = as_column_vector_or_scalar(y);
  using T_return = std::decay_t<decltype(f(x, y[0]))>;
  std::vector<T_return> result(y.size());
  Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
      = y_vec.unaryExpr([&f, &x](const auto& v) { return f(x, v); });
  return result;
}

/**
 * Specialisation for use with two nested containers (std::vectors).
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam Container1 Type of first std::vector to which functor is applied.
 * @tparam Container2 Type of second std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First std::vector input to which operation is applied.
 * @param y Second std::vector input to which operation is applied.
 * @param f functor to apply to std::vector inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <
    typename Container1, typename Container2, typename F,
    require_all_std_vector_vt<is_container_or_var_matrix, Container1, Container2>* = nullptr>
inline auto apply_scalar_binary(const Container1& x, const Container2& y, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  using T_return = plain_type_t<decltype(apply_scalar_binary(x[0], y[0], f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_binary(x[i], y[i], f);
  }
  return result;
}

/**
 * Specialisation for use when the first input is a nested std::vector and the
 * second is a scalar. The returned scalar type is deduced to allow for cases
 * where the input and return scalar types differ (e.g., functions implicitly
 * promoting integers).
 *
 * @tparam Container Type of std::vector to which functor is applied.
 * @tparam Scalar Type of scalar to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x std::vector input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename Container, typename Scalar, typename F,
          require_std_vector_vt<is_container_or_var_matrix, Container>* = nullptr,
          require_stan_scalar_t<Scalar>* = nullptr>
inline auto apply_scalar_binary(const Container& x, const Scalar& y, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_binary(x[0], y, f))>;
  size_t x_size = x.size();
  std::vector<T_return> result(x_size);
  for (size_t i = 0; i < x_size; ++i) {
    result[i] = apply_scalar_binary(x[i], y, f);
  }
  return result;
}

/**
 * Specialisation for use when the first input is a scalar and the second is a
 * nested std::vector. The returned scalar type is deduced to allow for cases
 * where the input and return scalar types differ (e.g., functions implicitly
 * promoting integers).
 *
 * @tparam Scalar Type of scalar to which functor is applied.
 * @tparam Scalar Type of std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Scalar input to which operation is applied.
 * @param y std::vector input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename Scalar, typename Container, typename F,
          require_stan_scalar_t<Scalar>* = nullptr,
          require_std_vector_vt<is_container_or_var_matrix, Container>* = nullptr>
inline auto apply_scalar_binary(const Scalar& x, const Container& y, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_binary(x, y[0], f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_binary(x, y[i], f);
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif
