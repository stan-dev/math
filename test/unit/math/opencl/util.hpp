#ifndef STAN_TEST_OPENCL_UTIL_HPP
#define STAN_TEST_OPENCL_UTIL_HPP

#include <stan/math.hpp>
#include <stan/math/opencl/rev/opencl.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>
#include <utility>
#include <string>
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
          require_all_not_st_var<T1, T2>* = nullptr>
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
template <typename T1, typename T2, require_all_rev_matrix_t<T1, T2>* = nullptr>
void expect_eq(const T1& a, const T2& b, const char* msg) {
  expect_eq(a.val(), b.val(), msg);
}
template <typename T>
void expect_eq(const std::vector<T>& a, const std::vector<T>& b,
               const char* msg) {
  EXPECT_EQ(a.size(), b.size());
  for (int i = 0; i < a.size(); i++) {
    expect_eq(a[i], b[i], msg);
  }
}
template <typename T1, typename T2,
          require_nonscalar_prim_or_rev_kernel_expression_t<T1>* = nullptr>
void expect_eq(const T1& a, const T2& b, const char* msg) {
  expect_eq(from_matrix_cl(a), b, msg);
}
template <typename T1, typename T2,
          require_nonscalar_prim_or_rev_kernel_expression_t<T2>* = nullptr>
void expect_eq(const T1& a, const T2& b, const char* msg) {
  expect_eq(a, from_matrix_cl(b), msg);
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

template <typename Functor>
void prim_rev_argument_combinations(Functor f) {
  f(std::make_tuple(), std::make_tuple());
}
template <typename Functor, typename Arg0, typename... Args>
void prim_rev_argument_combinations(const Functor& f, const Arg0& arg0,
                                    const Args&... args) {
  prim_rev_argument_combinations(
      [&f, &arg0](auto args_for_cpu, auto args_for_opencl) {
        constexpr size_t Size
            = std::tuple_size<std::decay_t<decltype(args_for_cpu)>>::value;
        return index_apply<Size>([&](auto... Is) {
          return f(
              std::forward_as_tuple(arg0, std::get<Is>(args_for_cpu)...),
              std::forward_as_tuple(arg0, std::get<Is>(args_for_opencl)...));
        });
      },
      args...);
  prim_rev_argument_combinations(
      [&](const auto& args_for_cpu, const auto& args_for_opencl) {
        constexpr size_t Size
            = std::tuple_size<std::decay_t<decltype(args_for_cpu)>>::value;
        return index_apply<Size>([&](auto... Is) {
          return f(std::make_tuple(var_argument(arg0),
                                   std::get<Is>(args_for_cpu)...),
                   std::make_tuple(var_argument(arg0),
                                   std::get<Is>(args_for_opencl)...));
        });
      },
      args...);
}

template <typename Functor, std::size_t... Is, typename... Args>
void compare_cpu_opencl_prim_rev_impl(const Functor& functor,
                                      std::index_sequence<Is...>,
                                      const Args&... args) {
  prim_rev_argument_combinations(
      [&functor](auto args_for_cpu, auto args_for_opencl) {
        auto res_cpu = functor(std::get<Is>(args_for_cpu)...);
        auto res_opencl
            = functor(opencl_argument(std::get<Is>(args_for_opencl))...);
        std::string signature = type_name<decltype(args_for_cpu)>().data();
        expect_eq(res_opencl, res_cpu,
                  ("CPU and OpenCL return values do not match for signature "
                   + signature + "!")
                      .c_str());
        var(recursive_sum(res_cpu) + recursive_sum(res_opencl)).grad();

        static_cast<void>(std::initializer_list<int>{
            (expect_adj_near(
                 std::get<Is>(args_for_opencl), std::get<Is>(args_for_cpu),
                 ("CPU and OpenCL adjoints do not match for argument "
                  + std::to_string(Is) + " for signature " + signature + "!")
                     .c_str()),
             0)...});

        set_zero_all_adjoints();
      },
      args...);
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
void compare_cpu_opencl_prim_rev(const Functor& functor, const Args&... args) {
  internal::compare_cpu_opencl_prim_rev_impl(
      functor, std::make_index_sequence<sizeof...(args)>{}, args...);
  recover_memory();
}

}  // namespace test
}  // namespace math
}  // namespace stan

#endif
