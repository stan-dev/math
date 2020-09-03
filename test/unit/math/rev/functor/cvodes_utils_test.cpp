#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <iostream>
#include <vector>

struct Inverse {
  template <typename T0, typename T_y>
  inline Eigen::Matrix<stan::return_type_t<T_y>, Eigen::Dynamic, 1> operator()(
      const T0& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
      std::ostream* msgs) const {
    Eigen::Matrix<T_y, Eigen::Dynamic, 1> out(1);
    out(0) = 1.0 / (y(0) - 1.0);
    return out;
  }
};

TEST(StanMath, cvodes_error_handler) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Ones(1);
  int t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  std::string msg
      = "CVODES: CVode Internal t = 0 and h = 0 are such that t + h = t on the "
        "next step";

  EXPECT_THROW_MSG(stan::math::ode_bdf(Inverse(), y0, t0, ts, nullptr),
                   std::domain_error, msg);
}
