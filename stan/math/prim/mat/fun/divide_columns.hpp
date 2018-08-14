#ifndef STAN_MATH_PRIM_MAT_FUN_DIVIDE_COLUMNS_HPP
#define STAN_MATH_PRIM_MAT_FUN_DIVIDE_COLUMNS_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/fun/divide.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Takes Stan data type vector[n] x[D] and divides each
 * dimension sequentially by each element of the vector
 *
 * @tparam T_x type of elements contained in vector x, scalar-type
 * @tparam T_v vector of elements
 * @tparam R   number of rows in the submatrix of x
 * @tparam C   number of columns in the submatrix of x
 *
 * @param x std::vector of elements
 * @param vec a vector type of elements
 * @throw std::invalid argument if D != length of vector
 *
 */
template <typename T_x, typename T_v, int R, int C>
inline typename std::vector<Eigen::Matrix<typename return_type<T_x, T_v,
                                                               double>::type,
                                          R, C>>
divide_columns(const std::vector<Eigen::Matrix<T_x, R, C>> &x,
               const std::vector<T_v> &vec) {
  size_t N = x.size();
  size_t D = x[0].size();
  check_size_match("divide_columns", "x dimension", D, "vector", vec.size());

  std::vector<Eigen::Matrix<
    typename return_type<T_x, T_v, double>::type, R, C>>
    out(N);
  for (size_t n = 0; n < N; ++n) {
    out[n].resize(D);
    for (size_t d = 0; d < D; ++d) {
      out[n][d] = divide(x[n][d], vec[d]);
    }
  }
  return out;
}
}  // namespace math
}  // namespace stan

#endif
