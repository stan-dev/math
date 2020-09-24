#ifndef STAN_TEST_OPENCL_UTIL_HPP
#define STAN_TEST_OPENCL_UTIL_HPP

#include <stan/math.hpp>
#include <stan/math/opencl/rev/opencl.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <gtest/gtest.h>
#include <utility>
#include <tuple>

namespace stan {
namespace math {
namespace test {

namespace internal {

template <typename T, require_stan_scalar_t<T>* = nullptr>
T opencl_argument(T&& x) {
  return x;
}

template <typename T, require_not_stan_scalar_t<T>* = nullptr>
auto opencl_argument(const T& x) {
  return to_matrix_cl(x);
}

template <typename T, require_st_same<T, int>* = nullptr>
T var_argument(T&& x) {
  return x;
}

template <typename T, require_not_st_same<T, int>* = nullptr>
auto var_argument(T&& x) {
  return to_var(x);
}

template <typename T, require_arithmetic_t<T>* = nullptr>
void expect_eq(T a, T b, const char* msg) {
  stan::test::expect_near_rel(msg, a, b);
}

void expect_eq(math::var a, math::var b, const char* msg) {
  stan::test::expect_near_rel(msg, a.val(), b.val());
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

template <typename T>
auto recursive_sum(const T& a) {
  return math::sum(a);
}

template <typename T>
auto recursive_sum(const std::vector<T>& a) {
  scalar_type_t<T> res = recursive_sum(a[0]);
  for (int i = 0; i < a.size(); i++) {
    res += recursive_sum(a[i]);
  }
  return res;
}

template <typename T, require_not_st_var<T>* = nullptr>
void expect_adj_near(const T& a, const T& b, const char* msg) {}

void expect_adj_near(var a, var b, const char* msg) {
  stan::test::expect_near_rel(msg, a.adj(), b.adj());
}

template <typename T1, typename T2, require_all_eigen_t<T1, T2>* = nullptr,
          require_vt_same<T1, T2>* = nullptr>
void expect_adj_near(const T1& a, const T2& b, const char* msg) {
  EXPECT_EQ(a.rows(), b.rows()) << msg;
  EXPECT_EQ(a.cols(), b.cols()) << msg;
  const auto& a_ref = math::to_ref(a);
  const auto& b_ref = math::to_ref(b);
  stan::test::expect_near_rel(msg, a_ref.adj(), b_ref.adj());
}

template <typename T>
void expect_adj_near(const std::vector<T>& a, const std::vector<T>& b,
                     const char* msg) {
  EXPECT_EQ(a.size(), b.size()) << msg;
  for (int i = 0; i < a.size(); i++) {
    expect_adj_near(a[i], b[i], msg);
  }
}

template <typename Functor, std::size_t... Is, typename... Args>
void compare_cpu_opencl_prim_rev_impl(Functor functor,
                                      std::index_sequence<Is...>,
                                      Args... args) {
  expect_eq(functor(opencl_argument(args)...), functor(args...),
            "CPU and OpenCL return values for prim do not match!");

  auto var_args_for_cpu = std::make_tuple(var_argument(args)...);
  auto var_args_for_opencl = std::make_tuple(var_argument(args)...);
  auto res_cpu = functor(std::get<Is>(var_args_for_cpu)...);
  auto res_opencl
      = functor(opencl_argument(std::get<Is>(var_args_for_opencl))...);
  expect_eq(res_opencl, res_cpu,
            "CPU and OpenCL return values for rev do not match!");

  (recursive_sum(res_cpu) + recursive_sum(res_opencl)).grad();

  static_cast<void>(std::initializer_list<int>{
      (expect_adj_near(
           std::get<Is>(var_args_for_opencl), std::get<Is>(var_args_for_cpu),
           "CPU and OpenCL adjoints do not match for argument " TO_STRING(
               Is) "!"),
       0)...});
}

}  // namespace internal

/**
 * Tests that given functor calculates same values and adjoints when given
 * arguments on CPU and OpenCL device.
 *
 * @tparam Functor type of the functor
 * @tparam Args types of the arguments
 * @param fucntor functor to test
 * @param args arguments to test the functor with. These should be just values
 * in CPU memory (no vars, no arguments on the OpenCL device).
 */
template <typename Functor, typename... Args>
void compare_cpu_opencl_prim_rev(Functor functor, Args... args) {
  internal::compare_cpu_opencl_prim_rev_impl(
      functor, std::make_index_sequence<sizeof...(args)>{}, args...);
}

}  // namespace test
}  // namespace math
}  // namespace stan

#endif
