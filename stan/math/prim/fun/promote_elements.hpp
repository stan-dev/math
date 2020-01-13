#ifndef STAN_MATH_PRIM_FUN_PROMOTE_ELEMENTS_HPP
#define STAN_MATH_PRIM_FUN_PROMOTE_ELEMENTS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {

/**
 * Struct with static function for elementwise type promotion.
 *
 * <p>This base implementation promotes one scalar value to another.
 *
 * @tparam T type of promoted element
 * @tparam S type of input element, must be assignable to T
 */
template <typename T, typename S>
struct promote_elements {
  /**
   * Return input element.
   *
   * @param u input of type S, assignable to type T
   * @returns input as type T
   */
  inline static T promote(const S& u) { return u; }
};

/**
 * Struct with static function for elementwise type promotion.
 *
 * <p>This specialization promotes scalar values of the same type.
 *
 * @tparam T type of elements
 */
template <typename T>
struct promote_elements<T, T> {
  /**
   * Return input element.
   *
   * @param u input of type T
   * @returns input as type T
   */
  inline static const T& promote(const T& u) { return u; }
};

/**
 * Struct with static function for elementwise type promotion.
 *
 * <p>This specialization promotes vector elements of different types
 * which must be compatible with promotion.
 *
 * @tparam T type of promoted elements
 * @tparam S type of input elements, must be assignable to T
 */
template <typename T, typename S>
struct promote_elements<std::vector<T>, std::vector<S> > {
  /**
   * Return input vector of type S as vector of type T.
   *
   * @param u vector of type S, assignable to type T
   * @returns vector of type T
   */
  inline static std::vector<T> promote(const std::vector<S>& u) {
    std::vector<T> t;
    t.reserve(u.size());
    for (size_t i = 0; i < u.size(); ++i) {
      t.push_back(promote_elements<T, S>::promote(u[i]));
    }
    return t;
  }
};

/**
 * Struct with static function for elementwise type promotion.
 *
 * <p>This specialization promotes vector elements of the same type.
 *
 * @tparam T type of elements
 */
template <typename T>
struct promote_elements<std::vector<T>, std::vector<T> > {
  /**
   * Return input vector.
   *
   * @param u vector of type T
   * @returns vector of type T
   */
  inline static const std::vector<T>& promote(const std::vector<T>& u) {
    return u;
  }
};

/**
 * Struct with static function for elementwise type promotion.
 *
 * <p>This specialization promotes matrix elements of different types
 * which must be compatible with promotion.
 *
 * @tparam T type of promoted elements
 * @tparam S type of input elements, must be assignable to T
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 */
template <typename T, typename S, int R, int C>
struct promote_elements<Eigen::Matrix<T, R, C>, Eigen::Matrix<S, R, C> > {
  /**
   * Return input matrix of type S as matrix of type T.
   *
   * @param u matrix of type S, assignable to type T
   * @returns matrix of type T
   */
  inline static Eigen::Matrix<T, R, C> promote(
      const Eigen::Matrix<S, R, C>& u) {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> t(u.rows(), u.cols());
    for (int i = 0; i < u.size(); ++i) {
      t(i) = promote_elements<T, S>::promote(u(i));
    }
    return t;
  }
};

/**
 * Struct with static function for elementwise type promotion.
 *
 * <p>This specialization promotes matrix elements of the same type.
 *
 * @tparam T type of elements in the matrices
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 */
template <typename T, int R, int C>
struct promote_elements<Eigen::Matrix<T, R, C>, Eigen::Matrix<T, R, C> > {
  /**
   * Return input matrix.
   *
   * @param u matrix of type T
   * @returns matrix of type T
   */
  inline static const Eigen::Matrix<T, R, C>& promote(
      const Eigen::Matrix<T, R, C>& u) {
    return u;
  }
};

}  // namespace math
}  // namespace stan

#endif
