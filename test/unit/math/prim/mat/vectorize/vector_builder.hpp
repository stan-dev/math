#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_VECTOR_BUILD_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_VECTOR_BUILD_HPP

#include <vector>

namespace test {
namespace math {

/**
 * Builder pattern for standard vectors.
 *
 * @tparam T Type of values held by builder.
 */
template <typename T>
struct vector_builder {
  std::vector<T> vec_;

  /**
   * Construct an empty vector builder.
   */
  vector_builder() {}

  /**
   * Add the specified element to the underlying vector
   * and return the builder for argument chaining.
   *
   * @param x Value to add.
   * @return Reference to this builder.
   */
  vector_builder& add(const T& x) {
    vec_.push_back(x);
    return *this;
  }

  /**
   * Return the underlying vector.
   *
   * @return Underlying vector.
   */
  std::vector<T> build() const { return vec_; }
};

}  // namespace math
}  // namespace test
#endif
