#ifndef STAN_MATH_PRIM_META_LENGTH_MVT_HPP
#define STAN_MATH_PRIM_META_LENGTH_MVT_HPP

#include <stanh/prim/meta/length_mvt.hpp>
#include <stanh/prim/fun/Eigen.hpp>
#include <stdexcept>
#include <vector>

namespace stan {

template <typename T, int R, int C>
size_t length_mvt(const Eigen::Matrix<T, R, C>& /* unused */) {
  return 1U;
}

template <typename T, int R, int C>
size_t length_mvt(const std::vector<Eigen::Matrix<T, R, C> >& x) {
  return x.size();
}

}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_META_LENGTH_MVT_HPP
#define STAN_MATH_PRIM_META_LENGTH_MVT_HPP

#include <stanh/prim/meta/length_mvt.hpp>
#include <stanh/prim/fun/Eigen.hpp>
#include <stdexcept>
#include <vector>

namespace stan {

template <typename T, int R, int C>
size_t length_mvt(const Eigen::Matrix<T, R, C>& /* unused */) {
  return 1U;
}

template <typename T, int R, int C>
size_t length_mvt(const std::vector<Eigen::Matrix<T, R, C> >& x) {
  return x.size();
}

}  // namespace stan
#endif
