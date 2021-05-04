#ifndef STAN_MATH_PRIM_FUN_DIVIDE_COLUMNS_HPP
#define STAN_MATH_PRIM_FUN_DIVIDE_COLUMNS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Takes Stan data type vector[n] x[D] and divides column
 * vector in x element-wise by the values in vec
 *
 * @tparam T_x Type of dividend
 * @tparam T_v Scalar type of divisor
 * @param x    std::vector of matrices
 * @param vec  std::vector of divisors
 * @throw std::invalid argument if D != length of vector
 */
template <typename T_x, typename T_v>
inline typename std::vector<
    Eigen::Matrix<return_type_t<T_x, T_v, double>, Eigen::Dynamic, 1>>
divide_columns(const std::vector<Eigen::Matrix<T_x, Eigen::Dynamic, 1>> &x,
               const std::vector<T_v> &vec) {
  const size_t N = x.size();
  const size_t D = x[0].size();
  check_size_match("divide_columns", "x dimension", D, "vector", vec.size());
  Eigen::Map<const Eigen::Array<T_v, Eigen::Dynamic, 1>> v_vec(&vec[0],
                                                               vec.size());

  std::vector<Eigen::Matrix<return_type_t<T_x, T_v, double>, Eigen::Dynamic, 1>>
      out(N);
  for (size_t n = 0; n < N; ++n) {
    out[n].resize(D);
    check_size_match("divide_columns", "x dimension", x[n].size(), "vector",
                     v_vec.size());
    out[n] = x[n].array() / v_vec.array();
  }
  return out;
}

}  // namespace math
}  // namespace stan

#endif
