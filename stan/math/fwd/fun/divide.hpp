#ifndef STAN_MATH_FWD_FUN_DIVIDE_HPP
#define STAN_MATH_FWD_FUN_DIVIDE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> divide(
    const Eigen::Matrix<fvar<T>, R, C>& v, const fvar<T>& c) {
  Eigen::Matrix<fvar<T>, R, C> res(v.rows(), v.cols());
  for (int i = 0; i < v.size(); i++)
    res(i) = v(i) / c;
  return res;
}

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> divide(
    const Eigen::Matrix<fvar<T>, R, C>& v, double c) {
  Eigen::Matrix<fvar<T>, R, C> res(v.rows(), v.cols());
  for (int i = 0; i < v.size(); i++)
    res(i) = v(i) / c;
  return res;
}

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> divide(const Eigen::Matrix<double, R, C>& v,
                                           const fvar<T>& c) {
  Eigen::Matrix<fvar<T>, R, C> res(v.rows(), v.cols());
  for (int i = 0; i < v.size(); i++)
    res(i) = v(i) / c;
  return res;
}

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> operator/(
    const Eigen::Matrix<fvar<T>, R, C>& v, const fvar<T>& c) {
  return divide(v, c);
}

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> operator/(
    const Eigen::Matrix<fvar<T>, R, C>& v, double c) {
  return divide(v, c);
}

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> operator/(
    const Eigen::Matrix<double, R, C>& v, const fvar<T>& c) {
  return divide(v, c);
}
}  // namespace math
}  // namespace stan
#endif
