#ifndef TEST_UNIT_MATH_MIX_MAT_SEQ_READER_HPP
#define TEST_UNIT_MATH_MIX_MAT_SEQ_READER_HPP

#include <stan/math/mix/mat.hpp>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {
namespace test {

/**
 * Structure for reading data from a serialized underlying
 * vector.
 *
 * @tparam scalar type of data
 */
template <typename T>
struct seq_reader {
  typedef Eigen::Matrix<T, -1, -1> matrix_t;
  typedef Eigen::Matrix<T, -1, 1> vector_t;
  typedef Eigen::Matrix<T, 1, -1> row_vector_t;

  /**
   * Underlying data.
   */
  vector_t data_;

  /**
   * Current position of read.
   */
  int pos_;

  /**
   * Construct a sequence reader from the specified data,
   * setting the position to the start of the sequence.
   *
   * @param data serialized data
   */
  explicit seq_reader(const vector_t& data) : data_(data), pos_(0) {}

  /**
   * Return the next scalar and advance the position.
   *
   * @return next data scalar
   * @throw std::logic_error if reading past end of reader
   */
  T next() {
    if (pos_ >= data_.size())
      throw std::logic_error("attempt to read past end of reader");
    return data_[pos_++];
  }

  /**
   * Read an object with size and container type determined by
   * the specified argument and return scalar type T determined
   * by the class instantiation.
   *
   * @param x value ignored, type specified to return a single
   * scalar
   * @return next scalar in the sequence
   */
  T read(double x) { return next(); }

  /**
   * Return a matrix read from the underlying reader whose
   * dimensions are given by the specified matrix argument.
   *
   * @tparam R static rows of matrix or vector
   * @tparam C static cols of matrix or vector
   * @param x argument specifying size and shape
   * @return matrix read from underlying stream of specified
   * size and shape
   */
  template <int R, int C>
  Eigen::Matrix<T, R, C> read(const Eigen::Matrix<double, R, C>& x) {
    Eigen::Matrix<T, R, C> y(x.rows(), x.cols());
    for (int i = 0; i < x.size(); ++i)
      y(i) = next();
    return y;
  }

  /**
   * Return a vector read from the underlying serialized reader
   * whose dimensions and types are given by specified argument.
   *
   * @tparam S argument array element type
   * @param x argument determing shape and size of result
   * @return return container of same shape and sizes as
   * argument read from the underlying scalar stream
   */
  template <typename S>
  typename promote_scalar_type<T, std::vector<S> >::type read(
      const std::vector<S>& x) {
    typename promote_scalar_type<T, std::vector<S> >::type y;
    y.reserve(x.size());
    for (size_t i = 0; i < x.size(); ++i)
      y.push_back(read(x[i]));
    return y;
  }
};

}  // namespace test
}  // namespace math
}  // namespace stan

#endif
