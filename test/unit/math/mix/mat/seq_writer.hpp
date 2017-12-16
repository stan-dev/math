#ifndef TEST_UNIT_MATH_MIX_MAT_SEQ_WRITER_HPP
#define TEST_UNIT_MATH_MIX_MAT_SEQ_WRITER_HPP

#include <stan/math/mix/mat.hpp>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {
namespace test {

/**
 * Structure holding a mutable vector of data into which
 * scalars, matrices, or arrays may be written in serial form.
 * The underlying data may be returned as an Eigen vector.  The
 * data written is all double based, but the accumulator may be
 * a type requiring promotion from double.
 *
 * @tparam T scalar type of data
 */
template <typename T>
struct seq_writer {
  /**
   * Accumulator for the data.
   */
  std::vector<T> data_;

  /**
   * Write the specified value into the data accumulator.
   *
   * @param x scalar to write
   */
  void write(double x) { data_.push_back(x); }

  /**
   * Write a serial version of the speicifed matrix into the
   * data accumulator in the order it is stored in the matrix.
   * Whether the argument is a matrix, vector or row vector is
   * determined by the template parameters.
   *
   * @tparam R static rows of matrix or vector
   * @tparam C static cols of matrix or vector
   * @param x matrix or vector
   */
  template <int R, int C>
  void write(const Eigen::Matrix<double, R, C>& x) {
    for (int i = 0; i < x.size(); ++i)
      write(x(i));
  }

  // TODO(carpenter): generalize to vector<T> and recurse

  /**
   * Write a vector of scalars into the underlying accumulator
   * in order.
   *
   * @param x vector to write
   */
  void write(const std::vector<double>& x) {
    for (size_t i = 0; i < x.size(); ++i)
      write(x[i]);
  }

  /**
   * Return a copy of the underlying accumulator as an Eigen
   * column vector.
   */
  Eigen::Matrix<T, -1, 1> vector() {
    Eigen::Matrix<T, -1, 1> y(data_.size());
    for (int i = 0; i < y.size(); ++i)
      y(i) = data_[i];
    return y;
  }
};

}  // namespace test
}  // namespace math
}  // namespace stan

#endif
