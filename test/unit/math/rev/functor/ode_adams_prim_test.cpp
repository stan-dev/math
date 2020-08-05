#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <iostream>
#include <vector>

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
  operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             std::ostream* msgs, const T_Args&... a) const {
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
  operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             std::ostream* msgs, const T2& a, const T3& b) const {
    Eigen::Matrix<stan::return_type_t<T1, T2, T3>, Eigen::Dynamic, 1> out(1);
    out << stan::math::cos((sum_(a) + sum_(b)) * t);
    return out;
  }
};

struct CosArgWrongSize {
  template <typename T0, typename T1, typename... T_Args>
  inline Eigen::Matrix<stan::return_type_t<T1, T_Args...>, Eigen::Dynamic, 1>
  operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             std::ostream* msgs, const T_Args&... a) const {
    std::vector<typename stan::return_type<T0, T_Args...>::type> vec
        = {sum_(a)...};
    Eigen::Matrix<stan::return_type_t<T1, T_Args...>, Eigen::Dynamic, 1> out(2);
    out << stan::math::cos(sum_(vec) * t), 0;
    return out;
  }
};

TEST(ode_adams_prim, y0_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  Eigen::VectorXd y0inf(1);
  Eigen::VectorXd y0NaN(1);
  Eigen::VectorXd y0_empty;
  y0NaN << stan::math::NOT_A_NUMBER;
  y0inf << stan::math::INFTY;
  int t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  EXPECT_NO_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0inf, t0, ts, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0NaN, t0, ts, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0_empty, t0, ts, nullptr, a),
               std::invalid_argument);
}

TEST(ode_adams_prim, t0_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  double t0inf = stan::math::INFTY;
  double t0NaN = stan::math::NOT_A_NUMBER;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  EXPECT_NO_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0inf, ts, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0NaN, ts, nullptr, a),
               std::domain_error);
}

TEST(ode_adams_prim, ts_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};
  std::vector<double> ts_repeat = {0.45, 0.45};
  std::vector<double> ts_lots = {0.45, 0.45, 1.1, 1.1, 2.0};
  std::vector<double> ts_empty = {};
  std::vector<double> ts_early = {-0.45, 0.2};
  std::vector<double> ts_decreasing = {0.45, 0.2};
  std::vector<double> tsinf = {stan::math::INFTY, 1.1};
  std::vector<double> tsNaN = {0.45, stan::math::NOT_A_NUMBER};

  double a = 1.5;

  std::vector<Eigen::VectorXd> out;
  EXPECT_NO_THROW(out
                  = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a));
  EXPECT_EQ(out.size(), ts.size());

  EXPECT_NO_THROW(
      out = stan::math::ode_adams(CosArg1(), y0, t0, ts_repeat, nullptr, a));
  EXPECT_EQ(out.size(), ts_repeat.size());
  EXPECT_MATRIX_FLOAT_EQ(out[0], out[1]);

  EXPECT_NO_THROW(
      out = stan::math::ode_adams(CosArg1(), y0, t0, ts_lots, nullptr, a));
  EXPECT_EQ(out.size(), ts_lots.size());

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts_empty, nullptr, a),
               std::invalid_argument);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts_early, nullptr, a),
               std::domain_error);

  EXPECT_THROW(
      stan::math::ode_adams(CosArg1(), y0, t0, ts_decreasing, nullptr, a),
      std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, tsinf, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, tsNaN, nullptr, a),
               std::domain_error);
}

TEST(ode_adams_prim, one_arg_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;
  double ainf = stan::math::INFTY;
  double aNaN = stan::math::NOT_A_NUMBER;

  std::vector<double> va = {a};
  std::vector<double> vainf = {ainf};
  std::vector<double> vaNaN = {aNaN};

  Eigen::VectorXd ea(1);
  ea << a;
  Eigen::VectorXd eainf(1);
  eainf << ainf;
  Eigen::VectorXd eaNaN(1);
  eaNaN << aNaN;

  std::vector<std::vector<double>> vva = {va};
  std::vector<std::vector<double>> vvainf = {vainf};
  std::vector<std::vector<double>> vvaNaN = {vaNaN};

  std::vector<Eigen::VectorXd> vea = {ea};
  std::vector<Eigen::VectorXd> veainf = {eainf};
  std::vector<Eigen::VectorXd> veaNaN = {eaNaN};

  EXPECT_NO_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, ainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, aNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, va));

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, vainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, vaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, ea));

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, eainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, eaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, vva));

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, vvainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, vvaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, vea));

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, veainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, veaNaN),
               std::domain_error);
}

TEST(ode_adams_prim, two_arg_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;
  double ainf = stan::math::INFTY;
  double aNaN = stan::math::NOT_A_NUMBER;

  std::vector<double> va = {a};
  std::vector<double> vainf = {ainf};
  std::vector<double> vaNaN = {aNaN};

  Eigen::VectorXd ea(1);
  ea << a;
  Eigen::VectorXd eainf(1);
  eainf << ainf;
  Eigen::VectorXd eaNaN(1);
  eaNaN << aNaN;

  std::vector<std::vector<double>> vva = {va};
  std::vector<std::vector<double>> vvainf = {vainf};
  std::vector<std::vector<double>> vvaNaN = {vaNaN};

  std::vector<Eigen::VectorXd> vea = {ea};
  std::vector<Eigen::VectorXd> veainf = {eainf};
  std::vector<Eigen::VectorXd> veaNaN = {eaNaN};

  EXPECT_NO_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, a));

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, ainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, aNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, va));

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, vainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, vaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, ea));

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, eainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, eaNaN),
               std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, vva));

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, vvainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, vvaNaN),
               std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, vea));

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, veainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a, veaNaN),
               std::domain_error);
}

TEST(ode_adams_prim, rhs_wrong_size_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  EXPECT_NO_THROW(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams(CosArgWrongSize(), y0, t0, ts, nullptr, a),
               std::invalid_argument);
}

TEST(ode_adams_prim, error_name) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double ainf = stan::math::INFTY;

  EXPECT_THROW_MSG(stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, ainf),
                   std::domain_error, "ode_adams");
}
