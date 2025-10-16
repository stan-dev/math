#ifndef STAN_TEST_OPENCL_UTIL_HPP
#define STAN_TEST_OPENCL_UTIL_HPP

#include <stan/math.hpp>
#include <stan/math/opencl/rev.hpp>
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
T opencl_argument(const T& x) {
  return x;
}
template <typename T, require_not_stan_scalar_t<T>* = nullptr>
auto opencl_argument(const T& x) {
  return to_matrix_cl(x);
}

template <typename T, require_t<std::is_integral<scalar_type_t<T>>>* = nullptr>
T var_argument(const T& x) {
  return x;
}
template <typename T,
          require_not_t<std::is_integral<scalar_type_t<T>>>* = nullptr,
          require_not_std_vector_t<T>* = nullptr>
auto var_argument(const T& x) {
  return to_var(x);
}
template <typename T,
          require_not_t<std::is_integral<scalar_type_t<T>>>* = nullptr,
          require_std_vector_t<T>* = nullptr>
auto var_argument(const T& x) {
  std::vector<decltype(var_argument(x[0]))> res;
  res.reserve(x.size());
  for (auto& i : x) {
    res.push_back(var_argument(i));
  }
  return res;
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
  expect_eq(from_matrix_cl<plain_type_t<T2>>(a), b, msg);
}
template <typename T1, typename T2,
          require_nonscalar_prim_or_rev_kernel_expression_t<T2>* = nullptr>
void expect_eq(const T1& a, const T2& b, const char* msg) {
  expect_eq(a, from_matrix_cl<plain_type_t<T1>>(b), msg);
}

template <typename T>
auto recursive_sum(const T& a) {
  return math::sum(a);
}
template <typename T>
auto recursive_sum(const std::vector<T>& a) {
  scalar_type_t<T> res = 0;
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
      [&f, &arg0](const auto& args1, const auto& args2) {
        constexpr size_t Size
            = std::tuple_size<std::decay_t<decltype(args1)>>::value;
        return index_apply<Size>([&](auto... Is) {
          return f(std::make_tuple(arg0, std::get<Is>(args1)...),
                   std::make_tuple(arg0, std::get<Is>(args2)...));
        });
      },
      args...);
  prim_rev_argument_combinations(
      [&](const auto& args1, const auto& args2) {
        constexpr size_t Size
            = std::tuple_size<std::decay_t<decltype(args1)>>::value;
        return index_apply<Size>([&](auto... Is) {
          return f(std::make_tuple(var_argument(arg0), std::get<Is>(args1)...),
                   std::make_tuple(var_argument(arg0), std::get<Is>(args2)...));
        });
      },
      args...);
}

template <typename Functor, std::size_t... Is, typename... Args>
void compare_cpu_opencl_prim_rev_impl(const Functor& functor,
                                      std::index_sequence<Is...>,
                                      const Args&... args) {
  prim_rev_argument_combinations(
      [&functor](const auto& args_for_cpu, const auto& args_for_opencl) {
        std::string signature = type_name<decltype(args_for_cpu)>().data();
        try {
          auto res_cpu = eval(functor(std::get<Is>(args_for_cpu)...));
          auto res_opencl = eval(
              functor(opencl_argument(std::get<Is>(args_for_opencl))...));
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
        } catch (...) {
          std::cerr << "exception thrown in signature " << signature << ":"
                    << std::endl;
          throw;
        }

        set_zero_all_adjoints();
      },
      args...);
}

template <typename FunctorCPU, typename FunctorCL, std::size_t... Is,
          typename... Args>
void compare_cpu_opencl_prim_rev_impl(const FunctorCPU& functorCPU,
                                      const FunctorCL& functorCL,
                                      std::index_sequence<Is...>,
                                      const Args&... args) {
  prim_rev_argument_combinations(
      [&functorCPU, &functorCL](const auto& args_for_cpu,
                                const auto& args_for_opencl) {
        auto res_cpu = eval(functorCPU(std::get<Is>(args_for_cpu)...));
        auto res_opencl = eval(
            functorCL(opencl_argument(std::get<Is>(args_for_opencl))...));
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

template <bool Condition, typename T, std::enable_if_t<Condition>* = nullptr>
auto to_vector_if(const T& x, std::size_t N) {
  return Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(N, x);
}
template <bool Condition, typename T, std::enable_if_t<!Condition>* = nullptr>
T to_vector_if(const T& x, std::size_t N) {
  return x;
}

using stan::math::rows;
template <typename T, require_not_container_t<T>* = nullptr>
int64_t rows(const T&) {
  return 1;
}
template <typename T, require_std_vector_t<T>* = nullptr>
int64_t rows(const T& x) {
  return x.size();
}

template <std::size_t I, typename Functor, std::size_t... Is, typename... Args>
void test_opencl_broadcasting_prim_rev_impl(const Functor& functor,
                                            std::index_sequence<Is...>,
                                            const Args&... args) {
  prim_rev_argument_combinations(
      [&functor, N = std::max({rows(args)...})](const auto& args_broadcast,
                                                const auto& args_vector) {
        auto res_scalar
            = eval(functor(opencl_argument(std::get<Is>(args_broadcast))...));
        std::string signature = type_name<decltype(args_broadcast)>().data();

        try {
          auto res_vec = eval(functor(opencl_argument(
              to_vector_if<Is == I>(std::get<Is>(args_vector), N))...));
          expect_eq(res_vec, res_scalar,
                    ("return values of broadcast and vector arguments do not "
                     "match for signature "
                     + signature + "!")
                        .c_str());
          try {
            var(recursive_sum(res_scalar) + recursive_sum(res_vec)).grad();
          } catch (...) {
            std::cerr << "throw in rev pass!" << std::endl;
            throw;
          }
        } catch (...) {
          std::cerr << "throw in signature: " << signature << "!" << std::endl;
          throw;
        }

        static_cast<void>(std::initializer_list<int>{
            (expect_adj_near(
                 std::get<Is>(args_vector), std::get<Is>(args_broadcast),
                 ("adjoints of broadcast and vector arguments do not match for "
                  "argument "
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
 * The functor must accept all possible combinations of converted `args`.
 * All std/row/col vectors and matrices are converted to `matrix_cl`. `double`
 * scalars can be converted to `var`s or left as `double`s. When converting
 * scalars to `double` in containers, `var_value<matrix_cl<double>>` is used
 * instead.
 *
 * @tparam Functor type of the functor
 * @tparam Args types of the arguments
 * @param fucntor functor to test
 * @param args arguments to test the functor with. These should be just values
 * in CPU memory (no vars, no arguments on the OpenCL device).
 */
template <typename Functor, typename... Args>
void compare_cpu_opencl_prim(const Functor& functor, const Args&... args) {
  auto res_cpu = eval(functor(args...));
  auto res_opencl = eval(functor(internal::opencl_argument(args)...));
  internal::expect_eq(res_cpu, res_opencl,
                      "CPU and OpenCL return values do not match!");
}

/**
 * Tests that given functor calculates same values and adjoints when given
 * arguments on CPU and OpenCL device.
 *
 * The functor must accept all possible combinations of converted `args`.
 * All std/row/col vectors and matrices are converted to `matrix_cl`. `double`
 * scalars can be converted to `var`s or left as `double`s. When converting
 * scalars to `double` in containers, `var_value<matrix_cl<double>>` is used
 * instead.
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

/**
 * Tests that given functor calculates same values and adjoints when given
 * arguments on CPU and OpenCL device.
 *
 * The functor must accept all possible combinations of converted `args`.
 * All std/row/col vectors and matrices are converted to `matrix_cl`. `double`
 * scalars can be converted to `var`s or left as `double`s. When converting
 * scalars to `double` in containers, `var_value<matrix_cl<double>>` is used
 * instead.
 *
 * @tparam FunctorCPU type of the CPU functor
 * @tparam FunctorCL type of the OpenCL functor
 * @tparam Args types of the arguments
 * @param fucntor functor to test
 * @param args arguments to test the functor with. These should be just values
 * in CPU memory (no vars, no arguments on the OpenCL device).
 */
template <typename FunctorCPU, typename FunctorCL, typename... Args>
void compare_cpu_opencl_prim_rev_separate(const FunctorCPU& functorCPU,
                                          const FunctorCL& fucntorCL,
                                          const Args&... args) {
  internal::compare_cpu_opencl_prim_rev_impl(
      functorCPU, fucntorCL, std::make_index_sequence<sizeof...(args)>{},
      args...);
  recover_memory();
}

/**
 * Tests that given functor can broadcast the `I`-th argument by comparing a
 * call with scalar and a call with column vector for that argument.
 *
 * The functor must accept all possible combinations of converted `args`.
 * All std/col vectors and matrices are converted to `matrix_cl`.
 * `double` scalars can be converted to `var`s or left as `double`s. When
 * converting scalars to `double` in containers, `var_value<matrix_cl<double>>`
 * is used instead.
 *
 * Warning: The scalar is broadcast to size equal to number of rows (or size in
 * case of `std::vector`) of the argument with the most rows
 *
 * @tparam Functor type of the functor
 * @tparam Args types of the arguments; `I`-th argument must be a scalar
 * @param fucntor functor to test
 * @param args arguments to test the functor with. These should be just values
 * in CPU memory (no vars, no arguments on the OpenCL device).
 */
template <std::size_t I, typename Functor, typename... Args,
          std::enable_if_t<(I < sizeof...(Args))>* = nullptr,
          require_stan_scalar_t<
              std::tuple_element_t<I, std::tuple<Args...>>>* = nullptr>
void test_opencl_broadcasting_prim_rev(const Functor& functor,
                                       const Args&... args) {
  internal::test_opencl_broadcasting_prim_rev_impl<I>(
      functor, std::make_index_sequence<sizeof...(args)>{}, args...);
  recover_memory();
}

}  // namespace test
}  // namespace math
}  // namespace stan

#endif
