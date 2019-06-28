#ifndef STANH_PRIM_METAAR_TYPE_HPP
#define STANH_PRIM_METAAR_TYPE_HPP
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <vector>
namespace stan {


/**
 * Metaprogram structure to determine the base scalar type
 * of a template argument.
 *
 * <p>This base class should be specialized for structured types.
 *
 * @tparam T Type of object.
 */
template <typename T>
struct scalar_type {
  typedef T type;
};

template <typename T>
struct scalar_type<T*> {
  typedef typename scalar_type<T>::type type;
};



/**
 * Template metaprogram defining the base scalar type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type of matrix.
 * @tparam R number of rows for matrix.
 * @tparam C number of columns for matrix.
 */
template <typename T, int R, int C>
struct scalar_type<Eigen::Matrix<T, R, C> > {
  typedef typename scalar_type<T>::type type;
};

/**
 * Template metaprogram defining the base scalar type of
 * values stored in a const Eigen matrix.
 *
 * @tparam T type of matrix.
 * @tparam R number of rows for matrix.
 * @tparam C number of columns for matrix.
 */
template <typename T, int R, int C>
struct scalar_type<const Eigen::Matrix<T, R, C> > {
  typedef typename scalar_type<T>::type type;
};

/**
 * Template metaprogram defining the base scalar type of
 * values stored in a referenced  Eigen matrix.
 *
 * @tparam T type of matrix.
 * @tparam R number of rows for matrix.
 * @tparam C number of columns for matrix.
 */
template <typename T, int R, int C>
struct scalar_type<Eigen::Matrix<T, R, C>&> {
  typedef typename scalar_type<T>::type type;
};

/**
 * Template metaprogram defining the base scalar type of
 * values stored in a referenced const Eigen matrix.
 *
 * @tparam T type of matrix.
 * @tparam R number of rows for matrix.
 * @tparam C number of columns for matrix.
 */
template <typename T, int R, int C>
struct scalar_type<const Eigen::Matrix<T, R, C>&> {
  typedef typename scalar_type<T>::type type;
};

/**
 * Template metaprogram defining the base scalar type of
 * values stored in an Eigen Block.
 *
 * @tparam T type of block.
 */
template <typename T>
struct scalar_type<Eigen::Block<T> > {
  typedef typename scalar_type<T>::type type;
};


template <typename T>
struct scalar_type<std::vector<T> > {
  typedef typename scalar_type<T>::type type;
};

template <typename T>
struct scalar_type<const std::vector<T> > {
  typedef typename scalar_type<T>::type type;
};

template <typename T>
struct scalar_type<std::vector<T>&> {
  typedef typename scalar_type<T>::type type;
};

template <typename T>
struct scalar_type<const std::vector<T>&> {
  typedef typename scalar_type<T>::type type;
};
#endif
