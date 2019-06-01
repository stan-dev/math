#ifndef TEST_UNIT_MATH_UTIL_HPP
#define TEST_UNIT_MATH_UTIL_HPP

#include <stan/math.hpp>
#include <string>
#include <vector>

namespace stan {
namespace test {

template <typename T>
Eigen::Matrix<T, -1, 1> to_eigen_vector(const std::vector<T>& x) {
  return Eigen::Map<const Eigen::Matrix<T, -1, 1>>(x.data(), x.size());
}

template <typename T, int R, int C>
std::vector<T> to_std_vector(const Eigen::Matrix<T, R, C>& x) {
  std::vector<T> y;
  y.reserve(x.size());
  for (int i = 0; i < x.size(); ++i)
    y.push_back(x(i));
  return y;
}

template <typename T>
struct deserializer {
  typedef T scalar_t;
  size_t position_;

  const std::vector<T> vals_;

  deserializer(const std::vector<T>& vals) : position_(0), vals_(vals) {}
  deserializer(const Eigen::Matrix<T, -1, 1>& v_vals)
      : position_(0), vals_(to_std_vector(v_vals)) {}

  template <typename U>
  T read(const U& x) {
    return vals_[position_++];
  }

  template <typename U>
  std::vector<T> read(const std::vector<U>& x) {
    std::vector<T> y;
    y.reserve(x.size());
    for (size_t i = 0; i < x.size(); ++i)
      y.push_back(read(x[i]));
    return y;
  }

  template <typename U, int R, int C>
  Eigen::Matrix<T, R, C> read(const Eigen::Matrix<U, R, C>& x) {
    Eigen::Matrix<T, R, C> y(x.rows(), x.cols());
    for (int i = 0; i < x.size(); ++i)
      y(i) = read(x(i));
    return y;
  }
};

template <typename T>
deserializer<T> to_deserializer(const std::vector<T>& vals) {
  return deserializer<T>(vals);
}
template <typename T>
deserializer<T> to_deserializer(const Eigen::Matrix<T, -1, 1>& vals) {
  return deserializer<T>(vals);
}

template <typename T>
struct serializer {
  typedef T scalar_t;
  std::vector<T> vals_;

  serializer() : vals_() {}

  // U assignable to T
  template <typename U>
  void write(const U& x) {
    vals_.push_back(x);
  }

  // U assignable to T
  template <typename U>
  void write(const std::vector<U>& x) {
    for (size_t i = 0; i < x.size(); ++i)
      write(x[i]);
  }

  // U assignable to T
  template <typename U, int R, int C>
  void write(const Eigen::Matrix<U, R, C>& x) {
    for (int i = 0; i < x.size(); ++i)
      write(x(i));
  }

  const std::vector<T>& array_vals() { return vals_; }
  const Eigen::Matrix<T, -1, 1>& vector_vals() {
    return to_eigen_vector(vals_);
  }
};

template <typename U>
void serialize_helper(serializer<U>& s) {}

template <typename U, typename T, typename... Ts>
void serialize_helper(serializer<U>& s, const T& x, const Ts... xs) {
  s.write(x);
  serialize_helper(s, xs...);
}

template <typename U, typename... Ts>
std::vector<U> serialize(const Ts... xs) {
  serializer<U> s;
  serialize_helper(s, xs...);
  return s.vals_;
}

template <typename T>
std::vector<typename scalar_type<T>::type> serialize_return(const T& x) {
  return serialize<typename scalar_type<T>::type>(x);
}

template <typename... Ts>
Eigen::VectorXd serialize_args(const Ts... xs) {
  return to_eigen_vector(serialize<double>(xs...));
}

}  // namespace test
}  // namespace stan

#endif
