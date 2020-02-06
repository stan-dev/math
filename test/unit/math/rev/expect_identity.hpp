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
  using stan::math::pack_params;

  Eigen::Matrix<double, sizeof...(xs), 1> x_vec(xs...);

  auto packed_f = pack_params<sizeof...(xs), F>(f);
  auto packed_g = pack_params<sizeof...(xs), G>(g);

  double value_f;
  Eigen::Matrix<double, Eigen::Dynamic, 1> grad_f;

  double value_g;
  Eigen::Matrix<double, Eigen::Dynamic, 1> grad_g;


  gradient(packed_f, x_vec, value_f, grad_f);
  gradient(packed_g, x_vec, value_g, grad_g);

  std::stringstream x_strstr;
  x_strstr << std::setprecision(24) << " for params: ";
  for(size_t i = 0; i < x_vec.size(); ++i) {
    if(i > 0) {
      x_strstr << ", ";
    }
    x_strstr << x_vec(i);
  }
  std::string x_str = x_strstr.str();

  std::stringstream msg_val;
  msg_val << msg << ": value" << x_str;
  expect_near_rel(msg_val.str(), value_f, value_g);
  for (size_t i = 0; i < grad_f.size(); ++i) {
    std::stringstream msg_grad;
    msg_grad << msg << ": grad(" << i << ")" << x_str;
    expect_near_rel(msg_grad.str(), grad_f(i), grad_g(i));
  }
}

}  // namespace test
}  // namespace stan

#endif
