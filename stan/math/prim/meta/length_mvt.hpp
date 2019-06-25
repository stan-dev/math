#ifndef STANH_PRIM_META_LENGTH_MVT_HPP
#define STANH_PRIM_META_LENGTH_MVT_HPP
#include <stdexcept>
#include <stan/math/prim/scal/meta/length_mvt.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stdexcept>
#include <vector>
namespace stan {




/**
 * length_mvt provides the length of a multivariate argument.
 *
 * This is the default template function. For any scalar type, this
 * will throw an std::invalid_argument exception since a scalar is not
 * a multivariate structure.
 *
 * @tparam T type to take length of. The default template function should
 *   only match scalars.
 * @throw std::invalid_argument since the type is a scalar.
 */
template <typename T>
size_t length_mvt(const T& /* unused */) {
  throw std::invalid_argument("length_mvt passed to an unrecognized type.");
  return 1U;
}




template <typename T, int R, int C>
size_t length_mvt(const Eigen::Matrix<T, R, C>& /* unused */) {
  return 1U;
}

template <typename T, int R, int C>
size_t length_mvt(const std::vector<Eigen::Matrix<T, R, C> >& x) {
  return x.size();
}

#endif
