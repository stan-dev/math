#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto lub_constrain_functor = [](const auto& a, const auto& b, const auto& c) {
  return stan::math::lub_constrain(a, b, c);
};
auto lub_constrain_functor2 = [](const auto& a, const auto& b, const auto& c) {
  using T_lp = stan::return_type_t<decltype(a), decltype(b), decltype(c)>;
  T_lp lp(4);
  if constexpr (!stan::is_constant<T_lp>::value) {
    stan::math::adjoint_of(lp) += 9;
  }
  return stan::math::lub_constrain(a, b, c, lp);
};
auto lub_constrain_functor3 = [](const auto& a, const auto& b, const auto& c) {
  using T_lp = stan::return_type_t<decltype(a), decltype(b), decltype(c)>;
  T_lp lp(4);
  stan::math::eval(stan::math::lub_constrain(a, b, c, lp));
  return lp;
};

TEST(OpenCLLubConstrain, prim_rev_values_small) {
  Eigen::VectorXd a(7);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4;
  Eigen::VectorXd b(7);
  b << -2.2, -0.8, 1, 4, -INFINITY, -INFINITY, 3;
  Eigen::VectorXd c(7);
  c << -2.1, 0.8, 4, 7, INFINITY, 2, INFINITY;
  double b_scal = -8;
  double c_scal = 8;

  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c_scal);
}

TEST(OpenCLLubConstrain, prim_rev_size_0) {
  int N = 0;

  Eigen::RowVectorXd a(N);
  Eigen::RowVectorXd b(N);
  Eigen::RowVectorXd c(N);
  double b_scal = -8;
  double c_scal = 8;

  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c_scal);
}

TEST(OpenCLLubConstrain, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, N);
  Eigen::MatrixXd c = b.array() + Eigen::ArrayXXd::Random(N, N) + 1.0;
  double b_scal = -1.5;
  double c_scal = 1.5;

  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c_scal);
}

TEST(OpenCLLubConstrain, rev_derivatives_tails) {
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  using stan::math::to_matrix_cl;
  using stan::math::var;
  using stan::math::var_value;
  // mpmath at 50 digits, lb = -2, ub = 3: dy/dx = 5 inv_logit(x) inv_logit(-x),
  // dy/dlb = inv_logit(-x), dlogJ/dx = inv_logit(-x) - inv_logit(x)
  Eigen::VectorXd x(7);
  x << -37, -20, 0, 20, 30, 37, 40;
  Eigen::VectorXd dx(7);
  dx << 4.2665238128720322e-16, 1.0305768069709247e-8, 1.25,
      1.0305768069709247e-8, 4.6788114844192117e-13, 4.2665238128720322e-16,
      2.1241771276457945e-17;
  Eigen::VectorXd dlb(7);
  dlb << 9.9999999999999991e-1, 9.9999999793884638e-1, 0.5,
      2.0611536181902036e-9, 9.357622968839299e-14, 8.5330476257440651e-17,
      4.248354255291589e-18;
  Eigen::VectorXd dlogj(7);
  dlogj << 9.9999999999999983e-1, 9.9999999587769276e-1, 0.0,
      -9.9999999587769276e-1, -9.9999999999981285e-1, -9.9999999999999983e-1,
      -9.9999999999999999e-1;
  Eigen::VectorXd lb = Eigen::VectorXd::Constant(7, -2.0);
  Eigen::VectorXd ub = Eigen::VectorXd::Constant(7, 3.0);
  for (bool use_lp : {false, true}) {
    var_value<matrix_cl<double>> x_cl = to_matrix_cl(x);
    var_value<matrix_cl<double>> lb_cl = to_matrix_cl(lb);
    var_value<matrix_cl<double>> ub_cl = to_matrix_cl(ub);
    var lp = 0;
    var s
        = use_lp
              ? stan::math::sum(
                  stan::math::lub_constrain(x_cl, lb_cl, ub_cl, lp))
              : stan::math::sum(stan::math::lub_constrain(x_cl, lb_cl, ub_cl));
    s.grad();
    Eigen::VectorXd x_adj = from_matrix_cl(x_cl.adj());
    Eigen::VectorXd lb_adj = from_matrix_cl(lb_cl.adj());
    for (int i = 0; i < x.size(); ++i) {
      EXPECT_NEAR(dx(i), x_adj(i), 1e-13 * dx(i)) << "dx at x = " << x(i);
      EXPECT_NEAR(dlb(i), lb_adj(i), 1e-13 * dlb(i)) << "dlb at x = " << x(i);
    }
    if (use_lp) {
      stan::math::set_zero_all_adjoints();
      lp.grad();
      x_adj = from_matrix_cl(x_cl.adj());
      for (int i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(dlogj(i), x_adj(i), 1e-13 * std::fabs(dlogj(i)))
            << "dlogJ/dx at x = " << x(i);
      }
    }
    stan::math::recover_memory();
  }
}

#endif
