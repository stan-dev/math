#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/rev/functor/idas_forward_system.hpp>
#include <stan/math/rev/functor/idas_integrator.hpp>

#include <nvector/nvector_serial.h>

#include <test/unit/util.hpp>
#include <gtest/gtest.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

struct chemical_kinetics {
  template <typename T0, typename TYY, typename TYP, typename TPAR>
  inline std::vector<stan::return_type_t<TYY, TYP, TPAR>> operator()(
      const T0& t_in, const std::vector<TYY>& yy, const std::vector<TYP>& yp,
      const std::vector<TPAR>& theta, const std::vector<double>& x_r,
      const std::vector<int>& x_i, std::ostream* msgs) const {
    if (yy.size() != 3 || yp.size() != 3)
      throw std::domain_error(
          "this function was called with inconsistent state");

    std::vector<stan::return_type_t<TYY, TYP, TPAR>> res(3);

    auto yy1 = yy.at(0);
    auto yy2 = yy.at(1);
    auto yy3 = yy.at(2);

    auto yp1 = yp.at(0);
    auto yp2 = yp.at(1);
    // auto yp3 = yp.at(2);

    auto p1 = theta.at(0);
    auto p2 = theta.at(1);
    auto p3 = theta.at(2);

    res[0] = yp1 + p1 * yy1 - p2 * yy2 * yy3;
    res[1] = yp2 - p1 * yy1 + p2 * yy2 * yy3 + p3 * yy2 * yy2;
    res[2] = yy1 + yy2 + yy3 - 1.0;

    // jacobian respect to yy

    return res;
  }
};

template <typename TYY, typename TYP, typename TPAR>
void idas_system_test(chemical_kinetics& f, std::vector<int> eq_id,
                      std::vector<TYY>& yy0, std::vector<TYP>& yp0,
                      std::vector<TPAR>& theta, N_Vector& res) {
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i(2, 0);
  std::ostream* msgs = 0;

  stan::math::idas_forward_system<chemical_kinetics, TYY, TYP, TPAR> dae{
      f, eq_id, yy0, yp0, theta, x_r, x_i, msgs};
  auto n = dae.n();
  EXPECT_EQ(n, yy0.size());
  EXPECT_EQ(dae.n_par(), theta.size());

  auto user_data = dae.to_user_data();
  auto yy = dae.nv_yy();
  auto yp = dae.nv_yp();
  for (size_t i = 0; i < n; ++i) {
    EXPECT_EQ(NV_Ith_S(yy, i), yy0[i]);
    EXPECT_EQ(NV_Ith_S(yp, i), yp0[i]);
  }

  auto residual = dae.residual();
  double t = 0.0;
  EXPECT_EQ(residual(t, yy, yp, res, user_data), 0);
  auto fval = f(t, yy0, yp0, theta, x_r, x_i, msgs);
  for (size_t i = 0; i < yy0.size(); ++i)
    EXPECT_EQ(NV_Ith_S(res, i), fval[i]);
}

template <typename TYY, typename TYP, typename TPAR>
void idas_forward_sen_test(chemical_kinetics& f, std::vector<int> eq_id,
                           std::vector<TYY>& yy0, std::vector<TYP>& yp0,
                           std::vector<TPAR>& theta, N_Vector& res,
                           N_Vector* ress) {
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i(2, 0);
  std::ostream* msgs = 0;

  stan::math::idas_forward_system<chemical_kinetics, TYY, TYP, TPAR> dae{
      f, eq_id, yy0, yp0, theta, x_r, x_i, msgs};

  auto yy = dae.nv_yy();
  auto yp = dae.nv_yp();
  auto yys = dae.nv_yys();
  auto yps = dae.nv_yps();
  auto n = dae.n();
  auto ns = dae.ns();

  auto user_data = dae.to_user_data();

  auto residual = dae.residual();
  double t = 0.0;
  EXPECT_EQ(residual(t, yy, yp, res, user_data), 0);

  if (dae.need_sens) {
    auto sensitivity_residual = dae.sensitivity_residual();
    auto temp1 = N_VNew_Serial(n);
    auto temp2 = N_VNew_Serial(n);
    auto temp3 = N_VNew_Serial(n);
    EXPECT_EQ(sensitivity_residual(ns, t, yy, yp, res, yys, yps, ress,
                                   user_data, temp1, temp2, temp3),
              0);
    N_VDestroy_Serial(temp1);
    N_VDestroy_Serial(temp2);
    N_VDestroy_Serial(temp3);
  }
}

TEST(IDAS_DAE_SYSTEM, idas_system_io) {
  chemical_kinetics f;
  std::vector<double> yy0{1.0, 0.0, 0.0};
  std::vector<double> yp0{0.1, 0.0, 0.0};
  std::vector<double> theta{0.040, 1.0e4, 3.0e7};
  const std::vector<int> eq_id{1, 1, 0};
  auto yy0_var = stan::math::to_var(yy0);
  auto yp0_var = stan::math::to_var(yp0);
  auto theta_var = stan::math::to_var(theta);

  {
    N_Vector res = N_VNew_Serial(yy0.size());
    idas_system_test(f, eq_id, yy0, yp0, theta, res);
    N_VDestroy_Serial(res);
  }

  {
    N_Vector res = N_VNew_Serial(yy0.size());
    idas_system_test(f, eq_id, yy0_var, yp0_var, theta_var, res);
    N_VDestroy_Serial(res);
  }

  {
    N_Vector res = N_VNew_Serial(yy0.size());
    idas_system_test(f, eq_id, yy0, yp0_var, theta_var, res);
    N_VDestroy_Serial(res);
  }
}

TEST(IDAS_DAE_SYSTEM, idas_forward_system_io) {
  chemical_kinetics f;
  std::vector<double> yy0{1.0, 0.0, 0.0};
  std::vector<double> yp0{0.1, 0.0, 0.0};
  std::vector<double> theta{0.040, 1.0e4, 3.0e7};
  const std::vector<int> eq_id{1, 1, 0};
  auto yy0_var = stan::math::to_var(yy0);
  auto yp0_var = stan::math::to_var(yp0);
  auto theta_var = stan::math::to_var(theta);
  const size_t n = yy0.size();

  {
    N_Vector res = N_VNew_Serial(n);
    N_Vector* ress = N_VCloneVectorArray(theta.size(), res);
    idas_forward_sen_test(f, eq_id, yy0, yp0, theta, res, ress);
    N_VDestroy_Serial(res);
    N_VDestroyVectorArray(ress, theta.size());
  }

  {
    N_Vector res = N_VNew_Serial(n);
    N_Vector* ress = N_VCloneVectorArray(theta.size(), res);
    idas_forward_sen_test(f, eq_id, yy0, yp0, theta_var, res, ress);
    EXPECT_EQ(NV_Ith_S(ress[0], 0), yy0[0]);
    EXPECT_EQ(NV_Ith_S(ress[0], 1), -yy0[0]);
    EXPECT_EQ(NV_Ith_S(ress[0], 2), 0);
    N_VDestroy_Serial(res);
    N_VDestroyVectorArray(ress, theta.size());
  }

  {
    size_t ns = n + n + theta.size();
    N_Vector res = N_VNew_Serial(n);
    N_Vector* ress = N_VCloneVectorArray(ns, res);
    idas_forward_sen_test(f, eq_id, yy0_var, yp0_var, theta_var, res, ress);
    EXPECT_EQ(NV_Ith_S(ress[n + n], 0), yy0[0]);
    EXPECT_EQ(NV_Ith_S(ress[n + n], 1), -yy0[0]);
    EXPECT_EQ(NV_Ith_S(ress[n + n], 2), 0);
    N_VDestroy_Serial(res);
    N_VDestroyVectorArray(ress, ns);
  }
}

TEST(IDAS_DAE_SYSTEM, idas_forward_system_general) {
  chemical_kinetics f;
  std::vector<double> yy0{0.5, 0.4, 0.2};
  std::vector<double> yp0{0.1, 0.2, 0.3};
  std::vector<double> theta{0.040, 1.0e1, 3.0e-1};
  const std::vector<int> eq_id{1, 1, 0};
  auto yy0_var = stan::math::to_var(yy0);
  auto yp0_var = stan::math::to_var(yp0);
  auto theta_var = stan::math::to_var(theta);
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i(2, 0);
  std::ostream* msgs = 0;

  stan::math::idas_forward_system<chemical_kinetics, stan::math::var,
                                  stan::math::var, stan::math::var>
      dae{f, eq_id, yy0_var, yp0_var, theta_var, x_r, x_i, msgs};
  size_t ns = dae.ns();
  size_t n = dae.n();
  auto yy = dae.nv_yy();
  auto yp = dae.nv_yp();
  auto yys = dae.nv_yys();
  auto yps = dae.nv_yps();

  auto p1 = theta[0];
  auto p2 = theta[1];
  auto p3 = theta[2];

  N_Vector res = N_VNew_Serial(n);
  N_Vector* ress = N_VCloneVectorArray(ns, res);
  for (size_t is = 0; is < ns; is++) {
    N_VConst(RCONST(1.0), yys[is]);
    N_VConst(RCONST(1.0), yps[is]);
  }
  auto user_data = dae.to_user_data();

  auto residual = dae.residual();
  double t = 0.0;
  EXPECT_EQ(residual(t, yy, yp, res, user_data), 0);
  auto fval = f(t, yy0, yp0, theta, x_r, x_i, msgs);
  for (size_t i = 0; i < n; ++i)
    EXPECT_EQ(NV_Ith_S(res, i), fval[i]);

  auto sensitivity_residual = dae.sensitivity_residual();
  auto temp1 = N_VNew_Serial(n);
  auto temp2 = N_VNew_Serial(n);
  auto temp3 = N_VNew_Serial(n);
  EXPECT_EQ(sensitivity_residual(ns, t, yy, yp, res, yys, yps, ress, user_data,
                                 temp1, temp2, temp3),
            0);
  Eigen::MatrixXd Jyy(n, n), Jyp(n, n), Jpar(n, ns);
  Jyy << p1, -p2 * yy0[2], -p2 * yy0[1], -p1, p2 * yy0[2] + 2 * p3 * yy0[1],
      p2 * yy0[1], 1, 1, 1;
  Jyp << 1, 0, 0, 0, 1, 0, 0, 0, 0;
  Jpar *= 0.0;
  Jpar.col(2 * n) << yy0[0], -yy0[0], 0;
  Jpar.col(2 * n + 1) << -yy0[1] * yy0[2], yy0[1] * yy0[2], 0;
  Jpar.col(2 * n + 2) << 0, yy0[1] * yy0[1], 0;
  Eigen::MatrixXd yys_mat(n, ns), yps_mat(n, ns);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < ns; ++j) {
      yys_mat(i, j) = 1.0;
      yps_mat(i, j) = 1.0;
    }
  }

  // check sensitivity_residual result
  auto r = Jyy * yys_mat + Jyp * yps_mat + Jpar;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < ns; ++j) {
      EXPECT_EQ(NV_Ith_S(ress[j], i), r(i, j));
    }
  }

  N_VDestroy_Serial(temp1);
  N_VDestroy_Serial(temp2);
  N_VDestroy_Serial(temp3);

  N_VDestroy_Serial(res);
  N_VDestroyVectorArray(ress, ns);
}

TEST(IDAS_DAE_SYSTEM, constructor_errors) {
  using stan::math::idas_forward_system;
  using stan::math::var;
  chemical_kinetics f;
  std::vector<double> yy0{0.5, 0.4, 0.2};
  std::vector<double> yp0{0.1, 0.2, 0.3};
  std::vector<double> theta{0.040, 1.0e1, 3.0e-1};
  const std::vector<int> eq_id{1, 1, 0};
  auto yy0_var = stan::math::to_var(yy0);
  auto yp0_var = stan::math::to_var(yp0);
  auto theta_var = stan::math::to_var(theta);
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i;
  std::ostream* msgs = 0;

  std::vector<double> bad_double{yy0};
  bad_double[0] = std::numeric_limits<double>::infinity();
  auto build_double = [&f, msgs, &x_i](const std::vector<int>& eq_id,
                                       const std::vector<double>& yy0,
                                       const std::vector<double>& yp0,
                                       const std::vector<double>& theta,
                                       const std::vector<double>& x_r) {
    idas_forward_system<chemical_kinetics, double, double, double> dae{
        f, eq_id, yy0, yp0, theta, x_r, x_i, msgs};
  };

  ASSERT_NO_THROW(build_double(eq_id, yy0, yp0, theta, x_r));

  EXPECT_THROW_MSG(build_double(eq_id, bad_double, yp0, theta, x_r),
                   std::domain_error, "initial state");

  EXPECT_THROW_MSG(build_double(eq_id, yp0, bad_double, theta, x_r),
                   std::domain_error, "derivative initial state");

  EXPECT_THROW_MSG(build_double(eq_id, yy0, yp0, bad_double, x_r),
                   std::domain_error, "parameter vector");

  EXPECT_THROW_MSG(build_double(eq_id, yy0, yp0, theta, bad_double),
                   std::domain_error, "continuous data");

  bad_double[0] = 0.0;
  bad_double.pop_back();
  EXPECT_THROW_MSG_WITH_COUNT(build_double(eq_id, bad_double, yp0, theta, x_r),
                              std::invalid_argument, "initial state", 2);

  std::vector<var> bad_var{std::numeric_limits<double>::infinity(), 1.0, 0.1};
  std::vector<var> empty_var;
  auto build_var = [&f, msgs, &x_i](const std::vector<int>& eq_id,
                                    const std::vector<var>& yy0,
                                    const std::vector<var>& yp0,
                                    const std::vector<var>& theta,
                                    const std::vector<double>& x_r) {
    idas_forward_system<chemical_kinetics, var, var, var> dae{
        f, eq_id, yy0, yp0, theta, x_r, x_i, msgs};
  };

  ASSERT_NO_THROW(build_var(eq_id, yy0_var, yp0_var, theta_var, x_r));

  EXPECT_THROW_MSG(build_var(eq_id, bad_var, yp0_var, theta_var, x_r),
                   std::domain_error, "initial state");

  EXPECT_THROW_MSG(build_var(eq_id, yp0_var, bad_var, theta_var, x_r),
                   std::domain_error, "derivative initial state");

  EXPECT_THROW_MSG(build_var(eq_id, yy0_var, yp0_var, bad_var, x_r),
                   std::domain_error, "parameter vector");

  EXPECT_THROW_MSG(build_var(eq_id, empty_var, yp0_var, theta_var, x_r),
                   std::invalid_argument, "initial state");

  EXPECT_THROW_MSG(build_var(eq_id, yp0_var, empty_var, theta_var, x_r),
                   std::invalid_argument, "derivative initial state");

  bad_var.erase(bad_var.begin());

  EXPECT_THROW_MSG(build_var(eq_id, yy0_var, bad_var, theta_var, x_r),
                   std::invalid_argument, "derivative initial state");

  std::vector<int> bad_eq_id{-1, 0, 0};
  EXPECT_THROW_MSG(build_var(bad_eq_id, yy0_var, yp0_var, theta_var, x_r),
                   std::domain_error, "derivative-algebra id");
  bad_eq_id = {1, 2, 0};
  EXPECT_THROW_MSG(build_var(bad_eq_id, yy0_var, yp0_var, theta_var, x_r),
                   std::domain_error, "derivative-algebra id");

  bad_eq_id.pop_back();
  EXPECT_THROW_MSG(build_var(bad_eq_id, yy0_var, yp0_var, theta_var, x_r),
                   std::invalid_argument, "derivative-algebra id");
}
