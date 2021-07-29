#ifndef TEST_UNIT_MATH_SERIALIZER_HPP
#define TEST_UNIT_MATH_SERIALIZER_HPP

#include <stan/math.hpp>
#include <complex>
#include <string>
#include <vector>

namespace stan {
namespace test {

/**
 * A class to store a sequence of values which can be deserialized
 * back into structured objects such as scalars, vectors, and
 * matrices.
 *
 * @tparam T type of scalars
 */
template <typename T>
struct deserializer {
  /**
   * Type of scalars in all objects.
   */
  typedef T scalar_t;

  /**
   * Current read position.
   */
  size_t position_;

  /**
   * The sequence of values to deserialize.
   */
  const std::vector<T> vals_;

  /**
   * Construct a deserializer from the specified sequence of values.
   *
   * @param vals values to deserialize
   */
  explicit deserializer(const std::vector<T>& vals)
      : position_(0), vals_(vals) {}

  /**
   * Construct a deserializer from the specified sequence of values.
   *
   * @param vals values to deserialize
   */
  explicit deserializer(const Eigen::Matrix<T, -1, 1>& v_vals)
      : position_(0), vals_(math::to_array_1d(v_vals)) {}

  /**
   * Read a scalar conforming to the shape of the specified argument,
   * here a scalar.  The specified argument is only used for its
   * shape---there is no relationship between the type of argument and
   * type of result.
   *
   * @tparam U type of pattern scalar
   * @param x pattern argument to determine result shape and size
   * @return deserialized value with shape and size matching argument
   */
  template <typename U, require_stan_scalar_t<U>* = nullptr,
            require_not_complex_t<U>* = nullptr>
  T read(const U& x) {
    return vals_[position_++];
  }

  /**
   * Read a complex number conforming to the shape of the specified
   * argument. The specified argument is only used for its
   * shape---there is no relationship between the value type of the
   * argument and the type of the result.
   *
   * @tparam U type of pattern value type
   * @param x pattern argument to determine result shape
   * @return deserialized value with shape and size matching argument
   */
  template <typename U>
  std::complex<T> read(const std::complex<U>& x) {
    T re = read(x.real());
    T im = read(x.imag());
    return {re, im};
  }

  /**
   * Read a standard vector conforming to the shape of the specified
   * argument, here a standard vector.  The specified argument is only
   * used for its shape---there is no relationship between the type of
   * argument and type of result.
   *
   * @tparam U type of pattern sequence elements
   * @param x pattern argument to determine result shape and size
   * @return deserialized value with shape and size matching argument
   */
  template <typename U, require_std_vector_t<U>* = nullptr,
            require_not_st_complex<U>* = nullptr>
  typename stan::math::promote_scalar_type<T, U>::type read(const U& x) {
    typename stan::math::promote_scalar_type<T, U>::type y;
    y.reserve(x.size());
    for (size_t i = 0; i < x.size(); ++i)
      y.push_back(read(x[i]));
    return y;
  }

  /**
   * Read a standard vector of std::complex variables conforming to the
   * shape of the specified argument, here a standard vector. The specified
   * argument is only used for its shape---there is no relationship between
   * the type of argument and type of result.
   *
   * @tparam U type of pattern sequence elements
   * @param x pattern argument to determine result shape and size
   * @return deserialized value with shape and size matching argument
   */
  template <typename U, require_std_vector_st<is_complex, U>* = nullptr>
  typename stan::math::promote_scalar_type<std::complex<T>, U>::type read(
      const U& x) {
    typename stan::math::promote_scalar_type<std::complex<T>, U>::type y;
    y.reserve(x.size());
    for (size_t i = 0; i < x.size(); ++i)
      y.push_back(read(x[i]));
    return y;
  }

  /**
   * Read a standard vector conforming to the shape of the specified
   * argument, here an Eigen matrix, vector, or row vector. The
   * specified argument is only used for its shape---there is no
   * relationship between the type of argument and type of result.
   *
   * @tparam U type of pattern scalar
   * @tparam R row specification for Eigen container
   * @tparam C column specification for Eigen container
   * @param x pattern argument to determine result shape and size
   * @return deserialized value with shape and size matching argument
   */
  template <typename U, int R, int C>
  Eigen::Matrix<T, R, C> read(const Eigen::Matrix<U, R, C>& x) {
    Eigen::Matrix<T, R, C> y(x.rows(), x.cols());
    for (int j = 0; j < x.cols(); ++j) {
      for (int i = 0; i < x.rows(); ++i) {
        y(i + j * x.rows()) = read(x(i, j));
      }
    }
    return y;
  }

  /**
   * Read a standard vector of std::complex variables conforming to the
   * shape of the specified argument, here an Eigen matrix, vector, or
   * row vector. The specified argument is only used for its shape---there
   * is no relationship between the type of argument and type of result.
   *
   * @tparam U type of pattern scalar
   * @tparam R row specification for Eigen container
   * @tparam C column specification for Eigen container
   * @param x pattern argument to determine result shape and size
   * @return deserialized value with shape and size matching argument
   */
  template <typename U, int R, int C>
  Eigen::Matrix<std::complex<T>, R, C> read(
      const Eigen::Matrix<std::complex<U>, R, C>& x) {
    Eigen::Matrix<std::complex<T>, R, C> y(x.rows(), x.cols());
    for (int j = 0; j < x.cols(); ++j) {
      for (int i = 0; i < x.rows(); ++i) {
        y(i + j * x.rows()) = read(x(i, j));
      }
    }
    return y;
  }
};

/**
 * A structure to serialize structures to an internal stored sequence
 * of scalars.
 *
 * @tparam T underlying scalar type
 */
template <typename T>
struct serializer {
  /**
   * Scalar type of serializer.
   */
  typedef T scalar_t;

  /**
   * Container for serialized values.
   */
  std::vector<T> vals_;

  /**
   * Construct a serializer.
   */
  serializer() : vals_() {}

  /**
   * Serialize the specified scalar.
   *
   * @tparam U type of specified scalar; must be assignable to T
   * @param x scalar to serialize
   */
  template <typename U, require_stan_scalar_t<U>* = nullptr>
  void write(const U& x) {
    vals_.push_back(x);
  }

  /**
   * Serialize the specified complex number.
   *
   * @tparam U value type of complex number; must be assignable to T
   * @param x complex number to serialize
   */
  template <typename U>
  void write(const std::complex<U>& x) {
    write(x.real());
    write(x.imag());
  }

  /**
   * Serialize the specified standard vector.
   *
   * @tparam U type of scalars; must be assignable to T
   * @param x vector to serialize
   */
  template <typename U>
  void write(const std::vector<U>& x) {
    for (size_t i = 0; i < x.size(); ++i)
      write(x[i]);
  }

  /**
   * Serialize the specified Eigen container.
   *
   * @tparam U type of scalars; must be assignable to T
   * @tparam R row specification of Eigen container
   * @tparam C column specification of Eigen container
   * @param x Eigen container to serialize.
   */
  template <typename U, require_eigen_t<U>* = nullptr>
  void write(const U& x) {
    for (int j = 0; j < x.cols(); ++j) {
      for (int i = 0; i < x.rows(); ++i) {
        write(x.coeff(i, j));
      }
    }
  }

  /**
   * Return the serialized values as a standard vector.
   *
   * @return serialized values
   */
  const std::vector<T>& array_vals() { return vals_; }

  /**
   * Return the serialized values as an Eigen vector.
   *
   * @return serialized values
   */
  const Eigen::Matrix<T, -1, 1>& vector_vals() {
    return math::to_vector(vals_);
  }
};

/**
 * Return a deserializer based on the specified values.
 *
 * @tparam T type of scalars in argument container and return
 * @param vals values to deserialize
 * @return deserializer based on specified values
 */
template <typename T>
deserializer<T> to_deserializer(const std::vector<T>& vals) {
  return deserializer<T>(vals);
}

/**
 * Return a deserializer based on the specified values.
 *
 * @tparam T type of scalars in argument container and return
 * @param vals values to deserialize
 * @return deserializer based on specified values
 */
template <typename T>
deserializer<T> to_deserializer(const Eigen::Matrix<T, -1, 1>& vals) {
  return deserializer<T>(vals);
}

template <typename T>
deserializer<T> to_deserializer(const std::complex<T>& vals) {
  return to_deserializer(std::vector<T>{vals.real(), vals.imag()});
}

template <typename U>
void serialize_helper(serializer<U>& s) {}

template <typename U, typename T, typename... Ts>
void serialize_helper(serializer<U>& s, const T& x, const Ts... xs) {
  s.write(x);
  serialize_helper(s, xs...);
}

/**
 * Serialize the specified sequence of objects, which all must have
 * scalar types assignable to the result scalar type.
 *
 * @tparam U type of scalar in result vector
 * @tparam Ts argument types
 * @param xs arguments
 * @return serialized form of arguments
 */
template <typename U, typename... Ts>
std::vector<U> serialize(const Ts... xs) {
  serializer<U> s;
  serialize_helper(s, xs...);
  return s.vals_;
}

/**
 * Serialized the specified single argument.
 *
 * @tparam T type of argument
 * @param x argument to serialize
 * @return serialized argument
 */
template <typename T>
std::vector<real_return_t<T>> serialize_return(const T& x) {
  return serialize<real_return_t<T>>(x);
}

/**
 * Serialize the specified sequence of structured objects with
 * double-based scalars into an Eigen vector of double values.
 *
 * @tparam Ts types of argument sequence
 * @param xs arguments to serialize
 * @return serialized form of arguments
 */
template <typename... Ts>
Eigen::VectorXd serialize_args(const Ts... xs) {
  return math::to_vector(serialize<double>(xs...));
}

}  // namespace test
}  // namespace stan

#endif
