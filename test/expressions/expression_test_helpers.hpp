#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <stan/math/fwd.hpp>
#include <vector>
#include <random>

namespace stan {
namespace test {

template <typename Scal>
struct counterOp {
  int* counter_;
  counterOp(int* counter) { counter_ = counter; }
  const Scal& operator()(const Scal& a) const {
    (*counter_)++;
    return a;
  }
};

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

template <typename T, require_integral_t<T>* = nullptr>
T make_arg(double value = 0.4, int size = 1) {
  return 1;
}
template <typename T, require_floating_point_t<T>* = nullptr>
T make_arg(double value = 0.4, int size = 1) {
  return value;
}
template <typename T, require_var_t<T>* = nullptr>
T make_arg(T value = 0.4, int size = 1) {
  return value;
}
template <typename T, require_fvar_t<T>* = nullptr>
T make_arg(T value, int size = 1) {
  return value;
}
template <typename T, require_fvar_t<T>* = nullptr>
T make_arg(double value = 0.4, int size = 1) {
  return {value, 0.5};
}
template <typename T, typename T_scalar = double,
          require_eigen_matrix_dynamic_t<T>* = nullptr>
T make_arg(T_scalar value = 0.4, int size = 1) {
  T res = T::Constant(size, size, make_arg<value_type_t<T>>(value));
  return res;
}
template <typename T, typename T_scalar = double,
          require_eigen_vector_t<T>* = nullptr>
T make_arg(T_scalar value = 0.4, int size = 1) {
  T res = T::Constant(size, make_arg<value_type_t<T>>(value));
  return res;
}
template <typename T, typename T_scalar = double,
          require_std_vector_t<T>* = nullptr>
T make_arg(T_scalar value = 0.4, int size = 1) {
  using V = value_type_t<T>;
  T res;
  for (int i = 0; i < size; i++) {
    res.push_back(make_arg<V>(value, size));
  }
  return res;
}
template <typename T, require_same_t<T, std::minstd_rand>* = nullptr>
T make_arg(double value = 0.4, int size = 1) {
  return std::minstd_rand(0);
}
template <typename T, require_same_t<T, std::ostream*>* = nullptr>
T make_arg(double value = 0.4, int size = 1) {
  return nullptr;
}

template <typename T>
T make_simplex(double value = 0.4, int size = 1) {
  return make_arg<T>(1.0 / size, size);
}

template <typename T>
T make_pos_definite_matrix(double value = 0.4, int size = 1) {
  T m1 = make_arg<T>(value, size);
  T m2 = m1 * stan::math::transpose(m1);
  m2.diagonal().array() += value;
  return m2;
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
  expect_eq(a.d(), b.d(), msg);
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
void expect_adj_eq(const T& a, const T& b, const char* msg = "expect_ad_eq") {}

void expect_adj_eq(math::var a, math::var b, const char* msg = "expect_ad_eq") {
  EXPECT_EQ(a.adj(), b.adj()) << msg;
}

template <typename T1, typename T2, require_all_eigen_t<T1, T2>* = nullptr,
          require_vt_same<T1, T2>* = nullptr>
void expect_adj_eq(const T1& a, const T2& b, const char* msg = "expect_ad_eq") {
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
                   const char* msg = "expect_ad_eq") {
  EXPECT_EQ(a.size(), b.size()) << msg;
  for (int i = 0; i < a.size(); i++) {
    expect_adj_eq(a[i], b[i], msg);
  }
}

void grad(stan::math::var& a) { a.grad(); }

#define TO_STRING_(x) #x
#define TO_STRING(x) TO_STRING_(x)
#define EXPECT_STAN_EQ(a, b) \
  stan::test::expect_eq(     \
      a, b, "Error in file: " __FILE__ ", on line: " TO_STRING(__LINE__))
#define EXPECT_STAN_ADJ_EQ(a, b) \
  stan::test::expect_adj_eq(     \
      a, b, "Error in file: " __FILE__ ", on line: " TO_STRING(__LINE__))

template <typename T, require_st_stan_scalar<T>* = nullptr>
scalar_type_t<T> sum_if_number(const T& a) {
  return recursive_sum(a);
}
template <typename T, require_not_st_stan_scalar<T>* = nullptr>
double sum_if_number(const T& a) {
  return 0;
}

struct test_functor {
  template <typename... T>
  auto operator()(T... args) const {
    using Ret_scal = return_type_t<T...>;
    return stan::test::make_arg<Eigen::Matrix<Ret_scal, Eigen::Dynamic, 1>>(
        math::sum(std::vector<Ret_scal>{
            static_cast<Ret_scal>(sum_if_number(args))...}));
  }
};

struct simple_eq_functor {
  template <typename T1, typename T2, typename T3, typename T4,
            require_eigen_vector_t<T1>* = nullptr,
            require_all_eigen_vector_t<T1, T2>* = nullptr>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T1& x, const T2& y, const T3& dat, const T4& dat_int,
             std::ostream* pstream__) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> z(1);
    z(0) = x(0) - y(0);
    return z;
  }
};

struct reduce_sum_test_functor {
  template <typename T0, typename T1>
  auto operator()(const T0& a, int, int, std::ostream*, const T1& b) {
    return math::sum(a) + math::sum(b);
  }
};

}  // namespace test

namespace math {

// for some functions signatures reported by stanc do not match ones in math. We
// use wrappers for such functions.
template <typename F, typename T_shared_param, typename T_job_param,
          typename T_x_r>
auto map_rect(const F& functor, const T_shared_param& shared_params,
              const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
                  job_params,
              const std::vector<std::vector<T_x_r>>& x_r,
              const std::vector<std::vector<int>>& x_i) {
  return map_rect<0, F>(shared_params, job_params, value_of(x_r), x_i);
}

template <typename F, typename T1, typename T2, typename T3, typename T4,
          typename T5>
auto algebra_solver_newton(const F& f, const T1& x, const T2& y,
                           const std::vector<double>& dat,
                           const std::vector<int>& dat_int,
                           T3 scaling_step_size, T4 function_tolerance,
                           T5 max_num_steps) {
  return algebra_solver_newton(
      f, x, y, dat, dat_int, nullptr, value_of(scaling_step_size),
      value_of(function_tolerance), value_of(max_num_steps));
}

template <typename F, typename T1, typename T2, typename T3, typename T4,
          typename T5>
auto algebra_solver(const F& f, const T1& x, const T2& y,
                    const std::vector<double>& dat,
                    const std::vector<int>& dat_int, T3 scaling_step_size,
                    T4 function_tolerance, T5 max_num_steps) {
  return algebra_solver(f, x, y, dat, dat_int, nullptr,
                        value_of(scaling_step_size),
                        value_of(function_tolerance), value_of(max_num_steps));
}

template <typename T1, typename T2>
auto reduce_sum(const T1& sliced, int grainsize, const T2& shared) {
  return reduce_sum<stan::test::reduce_sum_test_functor>(sliced, grainsize,
                                                         nullptr, shared);
}

// functions for testing the expression testing framework
template <typename T>
auto bad_no_expressions(const Eigen::Matrix<T, -1, -1>& a) {
  return a;
}

template <typename T>
auto bad_multiple_evaluations(const T& a) {
  return a + a;
}

/**
 * The bad_* functions are meant to be used when generating a failure
 */
template <typename T>
auto bad_wrong_value(const T& a) {
  if (std::is_same<T, plain_type_t<T>>::value) {
    return a(0, 0);
  }

  return a(0, 0) + 1;
}

template <typename T>
auto bad_wrong_derivatives(const T& a) {
  operands_and_partials<ref_type_t<T>> ops(a);
  if (!is_constant<T>::value && std::is_same<T, plain_type_t<T>>::value) {
    ops.edge1_.partials_[0] = 1234;
  }
  return ops.build(0);
}

}  // namespace math
}  // namespace stan
