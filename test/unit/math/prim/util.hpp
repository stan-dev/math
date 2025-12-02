#ifndef TEST_MATH_UNIT_PRIM_UTIL_HPP
#define TEST_MATH_UNIT_PRIM_UTIL_HPP

#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

namespace stan {
namespace test {
template <typename T>
struct always_true {
  static constexpr bool value = true;
};
template <typename T>
struct always_false {
  static constexpr bool value = false;
};

namespace internal {

template <typename T>
struct contains_floating_point {
  static constexpr bool value = std::is_floating_point_v<std::decay_t<T>>;
};
template <typename... Types>
struct contains_floating_point<std::tuple<Types...>> {
  static constexpr bool value
      = (contains_floating_point<std::decay_t<Types>>::value || ...);
};

template <typename T>
struct contains_std_vector {
  static constexpr bool value = stan::is_std_vector<std::decay_t<T>>::value;
};

template <typename... Types>
struct contains_std_vector<std::tuple<Types...>> {
  static constexpr bool value
      = (stan::is_std_vector<std::decay_t<Types>>::value || ...);
};
}  // namespace internal
template <typename... Types>
struct contains_floating_point {
  static constexpr bool value
      = (internal::contains_floating_point<std::decay_t<Types>>::value || ...);
};

template <typename... Types>
struct contains_std_vector {
  static constexpr bool value
      = (internal::contains_std_vector<std::decay_t<Types>>::value || ...);
};

template <Eigen::Index Idx, typename Tuple>
static constexpr bool is_const_ref_element_v = std::is_const_v<
    std::remove_reference_t<std::tuple_element_t<Idx, Tuple>>>&&
    std::is_reference_v<std::tuple_element_t<Idx, Tuple>>;

template <typename T>
struct is_fp_or_std_vector
    : std::bool_constant<
          internal::contains_floating_point<std::decay_t<T>>::value
          || internal::contains_std_vector<std::decay_t<T>>::value> {};

template <Eigen::Index Idx, typename Tuple>
static constexpr bool is_lvalue_ref_element_v
    = std::is_lvalue_reference_v<std::tuple_element_t<Idx, Tuple>>;

template <Eigen::Index Idx, typename Tuple>
static constexpr bool is_ref_element_v
    = std::is_reference_v<std::tuple_element_t<Idx, Tuple>>;

template <Eigen::Index Idx, typename Tuple>
static constexpr bool is_rvalue_element_v
    = std::is_rvalue_reference_v<std::tuple_element_t<Idx, Tuple>>;

template <Eigen::Index Idx, typename Tuple, typename T>
static constexpr bool is_same_tuple_element_v
    = std::is_same<std::decay_t<std::tuple_element_t<Idx, Tuple>>, T>::value;

namespace unit {

/**
 * Run a test that fails if the specified square matrix is not
 * symmetric.
 *
 * @param[in] a Matrix to test.
 */
inline void expect_symmetric(const Eigen::MatrixXd& a) {
  for (int j = 1; j < a.cols(); ++j)
    for (int i = 0; i < j; ++i)
      EXPECT_EQ(a(i, j), a(j, i)) << "failed symmetry at " << i << ", " << j;
}

/**
 * Return a randomly generated symmetric, positive-definite
 * matrix of the specified dimensionality using the specified
 * rng.
 *
 * @tparam RNG Class of random number generator.
 * @param[in] k Number of rows and columns in generated matrix.
 * @param[in, out] rng Random number generator.
 * @return Random k x k symmetric, positive-definite matrix.
 */
template <typename RNG>
Eigen::MatrixXd spd_rng(int k, RNG& rng) {
  using Eigen::MatrixXd;
  using stan::math::normal_rng;
  MatrixXd sigma = MatrixXd::Zero(k, k);
  for (int j = 0; j < k; ++j)
    for (int i = 0; i <= j; ++i)
      sigma(i, j) = normal_rng(0, 1, rng);
  for (int i = 0; i < k; ++i)
    sigma(i, i) *= sigma(i, i);               // pos. diagonal
  sigma = sigma.transpose() * sigma;          // reconstruct full matrix
  sigma = 0.5 * (sigma + sigma.transpose());  // symmetrize
  for (int i = 0; i < k; ++i)
    sigma(i, i) += 5;  // condition
  return sigma;
}
}  // namespace unit
}  // namespace test
}  // namespace stan
#endif
