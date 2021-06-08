#ifndef TEST_UNIT_MATH_PRIM_FUNCTOR_ODE_TEST_FUNCTORS_HPP
#define TEST_UNIT_MATH_PRIM_FUNCTOR_ODE_TEST_FUNCTORS_HPP

#include <stan/math/prim.hpp>

namespace stan {
namespace test {

template <typename T, stan::require_stan_scalar_t<T>* = nullptr>
T sum_(T arg) {
  return arg;
}

template <typename EigMat, stan::require_eigen_t<EigMat>* = nullptr>
auto sum_(EigMat&& arg) {
  return stan::math::sum(arg);
}

template <typename Vec, stan::require_std_vector_t<Vec>* = nullptr>
auto sum_(Vec&& arg) {
  stan::scalar_type_t<Vec> sum = 0;
  for (size_t i = 0; i < arg.size(); ++i) {
    sum += sum_(arg[i]);
  }
  return sum;
}

struct CosArg1 {
  template <typename T0, typename T1, typename... T_Args>
  inline Eigen::Matrix<stan::return_type_t<T1, T_Args...>, Eigen::Dynamic, 1>
  operator()(const T0& t, const T1& y, std::ostream* msgs,
             const T_Args&... a) const {
    std::vector<typename stan::return_type<T0, T_Args...>::type> vec
        = {sum_(a)...};
    Eigen::Matrix<stan::return_type_t<T1, T_Args...>, Eigen::Dynamic, 1> out(1);
    out << stan::math::cos(sum_(vec) * t);
    return out;
  }
};

struct Cos2Arg {
  template <typename T0, typename T1, typename T2, typename T3>
  inline Eigen::Matrix<stan::return_type_t<T1, T2, T3>, Eigen::Dynamic, 1>
  operator()(const T0& t, const T1& y, std::ostream* msgs, const T2& a,
             const T3& b) const {
    Eigen::Matrix<stan::return_type_t<T1, T2, T3>, Eigen::Dynamic, 1> out(1);
    out << stan::math::cos((sum_(a) + sum_(b)) * t);
    return out;
  }
};

struct CosArgWrongSize {
  template <typename T0, typename T1, typename... T_Args>
  inline Eigen::Matrix<stan::return_type_t<T1, T_Args...>, Eigen::Dynamic, 1>
  operator()(const T0& t, const T1& y, std::ostream* msgs,
             const T_Args&... a) const {
    std::vector<typename stan::return_type<T0, T_Args...>::type> vec
        = {sum_(a)...};
    Eigen::Matrix<stan::return_type_t<T1, T_Args...>, Eigen::Dynamic, 1> out(2);
    out << stan::math::cos(sum_(vec) * t), 0;
    return out;
  }
};

struct ayt {
  template <typename T0, typename T_y, typename T2>
  inline Eigen::Matrix<stan::return_type_t<T_y, T2>, Eigen::Dynamic, 1>
  operator()(const T0& t, const T_y& y, std::ostream* msgs, const T2& a) const {
    return -a * y * t;
  }
};

}  // namespace test
}  // namespace stan

#endif
