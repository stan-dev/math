#include <stan/math/rev/core.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

template <typename F, int R, int C>
void expect_ad2(const F& f, const Eigen::Matrix<double, R, C>& x) {
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::sum;
  
  Eigen::Matrix<var, R, C> xv = x;
  var_value<Eigen::Matrix<double, R, C>> vx = x;

  auto o = f(x);
  auto ov = f(xv);
  auto vo = f(vx);

  stan::math::set_zero_all_adjoints();
  sum(ov).grad();
  auto xv_adj = xv.adj().eval();
  stan::math::set_zero_all_adjoints();
  sum(vo).grad();
  auto vx_adj = vx.adj().eval();

  auto g = [&](const Eigen::VectorXd& flat_x) {
    Eigen::Matrix<double, R, C> x_shaped_variable(x.rows(), x.cols());
    for(size_t i = 0; i < flat_x.size(); ++i)
      x_shaped_variable(i) = flat_x(i);

    auto o = f(x_shaped_variable);

    return sum(o);
  };

  double tmp;
  Eigen::VectorXd flat_x(x.size());
  Eigen::VectorXd grad_x;
  for(size_t i = 0; i < flat_x.size(); ++i)
    flat_x(i) = x(i);
  stan::math::finite_diff_gradient_auto(g, flat_x, tmp, grad_x);
  Eigen::Matrix<double, R, C> x_shaped_grad(x.rows(), x.cols());
  for(size_t i = 0; i < grad_x.size(); ++i)
    x_shaped_grad(i) = grad_x(i);

  /*std::cout << x_shaped_grad << std::endl;
  std::cout << "-----" << std::endl;
  std::cout << vx_adj << std::endl;*/
   
  stan::test::expect_near_rel("matrix of vars function output", o, value_of(ov));
  stan::test::expect_near_rel("var matrix function output", o, value_of(vo));
  stan::test::expect_near_rel("matrix of vars argument adjoints", x_shaped_grad, xv_adj);
  stan::test::expect_near_rel("var matrix argument adjoints", x_shaped_grad, vx_adj);
}
