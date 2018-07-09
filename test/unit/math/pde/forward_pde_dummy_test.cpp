#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/arr/functor/forward_pde.hpp>
#include <stan/math/rev/arr/functor/forward_pde.hpp>
#include <gtest/gtest.h>
#include <sys/time.h>

#include <vector>
#include <iostream>
#include <iomanip>


/*
  test forward_pde func using hard-coded dummy map:
  y = f(theta) = { coef_ * theta[0] + coef_ * theta[1] ^ 2}
 */
class DummyPDEModel {
  const double coef_;
public:
  DummyPDEModel() : coef_(1.52) {}

  inline std::vector<std::vector<double> >
  operator()(const std::vector<double>& theta,
             const bool need_sens,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i,
             std::ostream* msgs = nullptr) const {
    std::vector<std::vector<double> > res;
    if(need_sens)
      res = {{coef_ * theta[0] + coef_ * theta[1] * theta[1],
              coef_,
              coef_ * 2.0 * theta[1]}};
    else
      res = {{coef_ * theta[0] + coef_ * theta[1] * theta[1]}};

    return res;
  }

  double coef() {return coef_;}
};

TEST(forward_pde, dummy_functor) {
  using stan::math::forward_pde;

  std::vector<double> x_r;
  std::vector<int> x_i;

  DummyPDEModel pde;
  std::vector<double> theta{1.0, 1.2};
  std::vector<double> qoi = forward_pde(pde, theta, x_r, x_i);
  ASSERT_EQ(qoi.size(), 1);
  ASSERT_FLOAT_EQ(qoi[0], pde.coef()*theta[0] + pde.coef()*theta[1]*theta[1]);

  const double p1 = 1.34;
  const double p2 = 2.81;
  std::vector<stan::math::var> theta_v{p1, p2};
  std::vector<stan::math::var> qoi_v =
    forward_pde(pde, theta_v, x_r, x_i);
  ASSERT_FLOAT_EQ(qoi_v[0].val(), pde.coef() * p1 + pde.coef() * p2 * p2);
  std::vector<double> g;
  stan::math::set_zero_all_adjoints();
  qoi_v[0].grad(theta_v, g);
  ASSERT_EQ(g.size(), theta.size());
  ASSERT_FLOAT_EQ(g[0], pde.coef());
  ASSERT_FLOAT_EQ(g[1], 2.0 * pde.coef() * p2);
}
