#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_BUILD_BINARY_VECTOR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_BUILD_BINARY_VECTOR_HPP

#include <vector>

template <typename F>
static inline std::vector<int> build_binary_vector1(std::vector<int> template_v,
                                                    int seed = -1) {
  return F::int_valid_inputs1();
}

template <typename F>
static inline std::vector<int> build_binary_vector2(std::vector<int> template_v,
                                                    int seed = -1) {
  return F::int_valid_inputs2();
}

template <typename F>
static inline std::vector<double> build_binary_vector1(
    std::vector<double> template_v, int seed = -1) {
  return F::valid_inputs1();
}

template <typename F>
static inline std::vector<double> build_binary_vector2(
    std::vector<double> template_v, int seed = -1) {
  return F::valid_inputs2();
}
#endif
