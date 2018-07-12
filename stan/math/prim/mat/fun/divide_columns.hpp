#ifndef STAN_MATH_PRIM_MAT_FUN_DIVIDE_COLUMNS_HPP
#define STAN_MATH_PRIM_MAT_FUN_DIVIDE_COLUMNS_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/fun/divide.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>

namespace stan {
namespace math {

/**
 * Takes Stan data type vector[n] x[D] and divides each
 * dimension by the specified scalar.
 *
 *
 * @tparam T_x type of elements contained in vector x, usually Eigen::Matrix
 * @tparam T_s type of element of scalar, usually double or var
 * @tparam R   number of rows in the submatrix of x
 * @tparam C   number of columns in the submatrix of x
 *
 * @param x std::vector of elements representing a matrix
 * @param scalar a scalar type of elementsx
 *
 */
template <typename T_x, typename T_s, int R, int C>
inline typename std::vector<Eigen::Matrix<T_x, R, C>> divide_columns(
    const std::vector<Eigen::Matrix<T_x, R, C>> &x, const T_s &scalar) {
  size_t N = x.size();
  std::vector<Eigen::Matrix<T_x, R, C>> out(N);
  for (size_t n = 0; n < N; ++n) {
    out[n] = divide(x[n], scalar);
  }
  return out;
}

/**
 * Takes Stan data type vector[n] x[D] and divides each
 * dimension sequentially by each element of the vector
 *
 * @tparam T_x type of elements contained in vector x, usually Eigen::Matrix
 * @tparam T_v vector of elements
 * @tparam R   number of rows in the submatrix of x
 * @tparam C   number of columns in the submatrix of x
 *
 * @param x std::vector of elements representing a matrix
 * @param scalar a scalar type of elementsx
 * @throw std::invalid argument if D != length of vector
 *
 */
template <typename T_x, typename T_v, int R, int C>
inline typename std::vector<Eigen::Matrix<T_x, R, C>> divide_columns(
    const std::vector<Eigen::Matrix<T_x, R, C>> &x,
    const std::vector<T_v> &vec) {
  size_t N = x.size();
  size_t D = x[0].size();
  check_size_match("divide_columns", "x dimension", D, "vector", vec.size());

  std::vector<Eigen::Matrix<T_x, R, C>> out(N);
  for (size_t n = 0; n < N; ++n) {
    for (size_t d = 0; d < D; ++d) {
      out[n].resize(D);
      out[n][d] = divide(x[n][d], vec[d]);
    }
  }
  return out;
}
}  // namespace math
}  // namespace stan

#endif
