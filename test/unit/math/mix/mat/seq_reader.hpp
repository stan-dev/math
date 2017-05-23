#ifndef TEST_UNIT_MATH_MIX_MAT_SEQ_READER_HPP
#define TEST_UNIT_MATH_MIX_MAT_SEQ_READER_HPP

#include <stan/math/mix/mat.hpp>
#include <stdexcept>

namespace stan {
  namespace math {
    namespace test {

      template <typename T>
      struct seq_reader {
        typedef Eigen::Matrix<T, -1, -1> matrix_t;
        typedef Eigen::Matrix<T, -1, 1> vector_t;
        typedef Eigen::Matrix<T, 1, -1> row_vector_t;

        vector_t data_;
        int pos_;

        seq_reader(const vector_t& data) : data_(data), pos_(0) { }

        T next() {
          if (pos_ >= data_.size())
            throw std::logic_error("attempt to read past end of reader");
          return data_[pos_++];
        }

        /**
         * Read an object with size and container type determined by
         * the specified argument and return scalar type T determined
         * by the class instantiation.
         */
        T read(double x) {
          return next();
        }

        template <int R, int C>
        matrix_t read(const Eigen::Matrix<double, R, C>& x) {
          Eigen::Matrix<T, R, C> y(x.rows(), x.cols());
          for (int i = 0; i < x.size(); ++i)
              y(i) = next();
          return y;
        }

        template <typename S>
        typename promote_scalar_type<T, std::vector<S> >::type
        read(const std::vector<S>& x) {
          typename promote_scalar_type<T, std::vector<S> >::type y;
          y.reserve(x.size());
          for (size_t i = 0; i < x.size(); ++i)
            y.push_back(read(x[i]));
          return y;
        }
      };

    }
  }
}

#endif
