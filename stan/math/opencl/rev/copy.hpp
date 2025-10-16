#ifndef STAN_MATH_OPENCL_REV_COPY_HPP
#define STAN_MATH_OPENCL_REV_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/vari.hpp>
#include <stan/math/opencl/rev/arena_type.hpp>
#include <stan/math/opencl/rev/to_arena.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>

#include <CL/opencl.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <type_traits>

namespace stan {
namespace math {

/** \ingroup opencl
 * Copies the source var containing Eigen matrices to destination var that has
 * data stored on the OpenCL device.
 *
 * @tparam T type of the Eigen matrix
 * @param a source Eigen matrix
 * @return var with a copy of the data on the OpenCL device
 */
template <typename T>
inline var_value<matrix_cl<value_type_t<T>>> to_matrix_cl(
    const var_value<T>& a) {
  return make_callback_var(to_matrix_cl(a.val()), [a](auto& res_vari) mutable {
    a.adj() += from_matrix_cl<plain_type_t<T>>(res_vari.adj());
  });
}

/** \ingroup opencl
 * Copies the source std::vector of vars to a destination var that has
 * data stored on the OpenCL device.
 *
 * @tparam T type of the std::vector
 * @param a source Eigen matrix
 * @return var with a copy of the data on the OpenCL device
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<value_type_t<T>>> to_matrix_cl(
    const std::vector<var_value<T>>& a) {
  return to_matrix_cl(
      Eigen::Map<const Eigen::Matrix<var_value<T>, Eigen::Dynamic, 1>>(
          a.data(), a.size()));
}

/** \ingroup opencl
 * Copies the source Eigen matrix of vars to
 * the destination matrix that is stored
 * on the OpenCL device.
 *
 * @tparam R Compile time rows of the Eigen matrix
 * @tparam C Compile time columns of the Eigen matrix
 * @param src source Eigen matrix
 * @return matrix_cl with a copy of the data in the source matrix
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline var_value<matrix_cl<value_type_t<value_type_t<T>>>> to_matrix_cl(
    const T& src) {
  arena_t<T> src_stacked = src;

  return make_callback_var(
      to_matrix_cl(src_stacked.val()), [src_stacked](auto& res_vari) mutable {
        src_stacked.adj() += from_matrix_cl<
            Eigen::Matrix<double, T::RowsAtCompileTime, T::ColsAtCompileTime>>(
            res_vari.adj());
      });
}

/** \ingroup opencl
 * Copies the source vector of Eigen matrices of vars to
 * the destination matrix that is stored
 * on the OpenCL device. Each element of the vector is stored into one column of
 * the returned matrix_cl.
 *
 * @param src source vector of Eigen matrices
 * @return matrix_cl with a copy of the data in the source matrix
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline var_value<matrix_cl<value_type_t<value_type_t<T>>>> to_matrix_cl(
    const std::vector<T>& src) {
  auto src_stacked = to_arena(src);

  return make_callback_var(
      to_matrix_cl(value_of(src_stacked)),
      [src_stacked](auto& res_vari) mutable {
        Eigen::MatrixXd adj = from_matrix_cl(res_vari.adj());
        for (int i = 0; i < src_stacked.size(); i++) {
          src_stacked[i].adj()
              += Eigen::Map<plain_type_t<decltype(src_stacked[i].adj())>>(
                  adj.data() + adj.rows() * i, src_stacked[i].rows(),
                  src_stacked[i].cols());
        }
      });
}

/** \ingroup opencl
 * Copies the source var that has data stored on the OpenCL device to
 * destination var containing Eigen matrix.
 *
 * @tparam T_dst destination type
 * @tparam T type of the matrix or expression on the OpenCL device
 * @param a source matrix_cl or expression
 * @return var with a copy of the data on the host
 */
template <typename T_dst, typename T,
          require_var_vt<is_eigen, T_dst>* = nullptr,
          require_all_kernel_expressions_t<T>* = nullptr>
inline T_dst from_matrix_cl(const var_value<T>& a) {
  return make_callback_var(
      from_matrix_cl<Eigen::Matrix<double, T_dst::RowsAtCompileTime,
                                   T_dst::ColsAtCompileTime>>(a.val()),
      [a](auto& res_vari) mutable { a.adj() += to_matrix_cl(res_vari.adj()); });
}

/** \ingroup opencl
 * Copies the source var that has data stored on the OpenCL device to
 * destination Eigen matrix containing vars.
 *
 * @tparam T_dst destination type
 * @tparam T type of the matrix or expression on the OpenCL device
 * @param a source matrix_cl or expression
 * @return var with a copy of the data on the host
 */
template <typename T_dst, typename T,
          require_eigen_vt<is_var, T_dst>* = nullptr,
          require_all_kernel_expressions_t<T>* = nullptr>
inline T_dst from_matrix_cl(const var_value<T>& a) {
  arena_t<T_dst> res
      = from_matrix_cl<Eigen::Matrix<double, T_dst::RowsAtCompileTime,
                                     T_dst::ColsAtCompileTime>>(a.val());
  reverse_pass_callback(
      [a, res]() mutable { a.adj() += to_matrix_cl(res.adj()); });
  return res;
}

/** \ingroup opencl
 * Copies the source var that has data stored on the OpenCL device to
 * destination `std::vector` containing vars.
 *
 * @tparam T_dst destination type
 * @tparam T type of the matrix or expression on the OpenCL device
 * @param a source matrix_cl or expression
 * @return var with a copy of the data on the host
 */
template <typename T_dst, typename T,
          require_std_vector_vt<is_var, T_dst>* = nullptr,
          require_all_stan_scalar_t<value_type_t<T_dst>>* = nullptr,
          require_all_kernel_expressions_t<T>* = nullptr>
inline T_dst from_matrix_cl(const var_value<T>& a) {
  check_size_match("from_matrix_cl<std::vector<var>>", "src.cols()", a.cols(),
                   "dst.cols()", 1);
  std::vector<double> val = from_matrix_cl<std::vector<double>>(a.val());
  arena_t<T_dst> res(val.begin(), val.end());
  reverse_pass_callback([a, res]() mutable {
    a.adj() += to_matrix_cl(as_column_vector_or_scalar(res).adj());
  });
  return {res.begin(), res.end()};
}

/** \ingroup opencl
 * Copies the source var that has data stored on the OpenCL device to
 * destination std::vector containing either Eigen vectors of vars or vars
 * containing Eigen vectors.
 *
 * @tparam T_dst destination type
 * @tparam T type of the matrix or expression on the OpenCL device
 * @param a source matrix_cl or expression
 * @return var with a copy of the data on the host
 */
template <typename T_dst, typename T, require_std_vector_t<T_dst>* = nullptr,
          require_rev_vector_t<value_type_t<T_dst>>* = nullptr,
          require_all_kernel_expressions_t<T>* = nullptr>
inline T_dst from_matrix_cl(const var_value<T>& a) {
  Eigen::MatrixXd val = from_matrix_cl(a.val());
  arena_t<T_dst> res;
  res.reserve(a.cols());
  for (int i = 0; i < a.cols(); i++) {
    res.emplace_back(val.col(i));
  }
  reverse_pass_callback([a, res]() mutable {
    Eigen::MatrixXd adj(a.rows(), a.cols());
    for (int i = 0; i < a.cols(); i++) {
      adj.col(i) = res[i].adj();
    }
    a.adj() += to_matrix_cl(adj);
  });
  return {res.begin(), res.end()};
}

/** \ingroup opencl
 * Copies the source var that has data stored on the OpenCL device to
 * destination Eigen matrix containing vars.
 *
 * @tparam T type of the matrix or expression on the OpenCL device
 * @param src source matrix_cl or expression
 * @return var with a copy of the data on the host
 */
template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
auto from_matrix_cl(const var_value<T>& src) {
  return from_matrix_cl<var_value<Eigen::MatrixXd>>(src);
}

}  // namespace math
}  // namespace stan
#endif
#endif
