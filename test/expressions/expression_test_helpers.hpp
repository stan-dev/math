#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <stan/math/fwd.hpp>
#include <vector>

namespace stan {
namespace test {

template<typename T, require_integral_t<T>* = nullptr> T make_arg() { return 1; }
template<typename T, require_floating_point_t<T>* = nullptr> T make_arg() { return 0.4; }
template<typename T, require_autodiff_t<T>* = nullptr> T make_arg() { return 0.4; }
template<typename T, require_eigen_vector_t<T>* = nullptr> T make_arg() {
  T res(3);
  res << 0.1, 0.2, 0.3;
  return res;
}
template<typename T, require_eigen_matrix_t<T>* = nullptr> T make_arg() {
  T res(3, 3);
  res << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9;
  return res;
}
template<typename T, require_std_vector_t<T>* = nullptr> T make_arg() {
  using V = value_type_t<T>;
  V tmp = make_arg<V>();
  T res;
  res.push_back(tmp);
  res.push_back(tmp);
  res.push_back(tmp);
  return res;
}

template <typename T, require_arithmetic_t<T>* = nullptr>
void expect_eq(T a, T b, const char* msg) {
  EXPECT_EQ(a, b) << msg;
}

void expect_eq(math::var a, math::var b, const char* msg) {
  EXPECT_EQ(a.val(), b.val()) << msg;
}

template <typename T, require_arithmetic_t<T>* = nullptr>
void expect_eq(math::fvar<T> a, math::fvar<T> b, const char* msg) {
  expect_eq(a.val(), b.val(), msg);
  expect_eq(a.tangent(), b.tangent(), msg);
}

template <typename T1, typename T2, require_all_eigen_t<T1, T2>* = nullptr,
          require_vt_same<T1, T2>* = nullptr>
void expect_eq(const T1& a, const T2& b, const char* msg) {
  EXPECT_EQ(a.rows(), b.rows()) << msg;
  EXPECT_EQ(a.cols(), b.cols()) << msg;
  const auto& a_ref = math::to_ref(a);
  const auto& b_ref = math::to_ref(b);
  for (int i = 0; i < a.rows(); i++) {
    for (int j = 0; j < a.cols(); j++) {
      expect_eq(a_ref(i, j), b_ref(i, j), msg);
    }
  }
}

template <typename T>
void expect_eq(const std::vector<T>& a, const std::vector<T>& b,
               const char* msg) {
  EXPECT_EQ(a.size(), b.size());
  for (int i = 0; i < a.size(); i++) {
    expect_eq(a[i], b[i], msg);
  }
}

template <typename T, require_not_st_var<T>* = nullptr>
void expect_adj_eq(const T& a, const T& b, const char* msg) {}

void expect_adj_eq(math::var a, math::var b, const char* msg) {
  EXPECT_EQ(a.adj(), b.adj()) << msg;
}

template <typename T1, typename T2, require_all_eigen_t<T1, T2>* = nullptr,
          require_vt_same<T1, T2>* = nullptr>
void expect_adj_eq(const T1& a, const T2& b, const char* msg) {
  EXPECT_EQ(a.rows(), b.rows()) << msg;
  EXPECT_EQ(a.cols(), b.cols()) << msg;
  const auto& a_ref = math::to_ref(a);
  const auto& b_ref = math::to_ref(b);
  for (int i = 0; i < a.rows(); i++) {
    for (int j = 0; j < a.cols(); j++) {
      expect_adj_eq(a_ref(i, j), b_ref(i, j), msg);
    }
  }
}

template <typename T>
void expect_adj_eq(const std::vector<T>& a, const std::vector<T>& b,
                     const char* msg) {
  EXPECT_EQ(a.size(), b.size()) << msg;
  for (int i = 0; i < a.size(); i++) {
    expect_adj_eq(a[i], b[i], msg);
  }
}

#define TO_STRING_(x) #x
#define TO_STRING(x) TO_STRING_(x)
#define EXPECT_STAN_EQ(a, b) \
  stan::test::expect_eq(      \
      a, b, "Error in file: " __FILE__ ", on line: " TO_STRING(__LINE__))
#define EXPECT_STAN_ADJ_EQ(a, b) \
  stan::test::expect_adj_eq(      \
      a, b, "Error in file: " __FILE__ ", on line: " TO_STRING(__LINE__))

}  // namespace math
}  // namespace stan
