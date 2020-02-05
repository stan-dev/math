#ifndef TEST_UNIT_MATH_REV_EXPECT_IDENTITY_HPP
#define TEST_UNIT_MATH_REV_EXPECT_IDENTITY_HPP

#include <stan/math/rev.hpp>
#include <test/unit/math/relative_tolerance.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <string>

namespace stan {
namespace test {

template <typename F, typename G, typename... Ts>
void expect_identity(const std::string& msg, const F& f, const G& g, Ts... xs) {
  using stan::math::gradient;

  double value_f;
  Eigen::Matrix<double, Eigen::Dynamic, 1> grad_f;

  double value_g;
  Eigen::Matrix<double, Eigen::Dynamic, 1> grad_g;

  gradient(f, value_f, grad_f, xs...);
  gradient(g, value_g, grad_g, xs...);

  std::stringstream msg_val;
  msg_val << msg << ": value";
  expect_near_rel(msg_val.str(), value_f, value_g);  
  for(size_t i = 0; i < grad_f.size(); ++i) {
    std::stringstream msg_grad;
    msg_grad << msg << ": grad(" << i << ")";
    expect_near_rel(msg_grad.str(), grad_f(i), grad_g(i));
  }

}

}
}

#endif