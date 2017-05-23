#ifndef TEST_UNIT_MATH_MIX_MAT_SEQ_WRITER_HPP
#define TEST_UNIT_MATH_MIX_MAT_SEQ_WRITER_HPP

#include <stan/math/mix/mat.hpp>
#include <stdexcept>

namespace stan {
  namespace math {
    namespace test {

      template <typename T>
      struct seq_writer {
        std::vector<T> data_;



        void write(double x) {
          data_.push_back(x);
        }

        template <int R, int C>
        void write(const Eigen::Matrix<double, R, C>& x) {
          for (int i = 0; i < x.size(); ++i)
            write(x(i));
        }

        void write(const std::vector<double>& x) {
          for (size_t i = 0; i < x.size(); ++i)
            write(x[i]);
        }

        Eigen::Matrix<T, -1, 1> vector() {
          Eigen::Matrix<T, -1, 1> y(data_.size());
          for (int i = 0; i < y.size(); ++i)
            y(i) = data_[i];
          return y;
        }
      };

    }
  }
}

#endif
