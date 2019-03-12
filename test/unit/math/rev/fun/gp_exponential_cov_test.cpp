#include <stan/math/rev.hpp>
#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>
#include <string>
#include <vector>

template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
std::string pull_msg(std::vector<T_x1> x1, std::vector<T_x2> x2, T_sigma sigma,
                     T_l l) {
  std::string message;
  try {
    stan::math::gp_exponential_cov(x1, x2, sigma, l);
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

template <typename T_x1, typename T_sigma, typename T_l>
std::string pull_msg(std::vector<T_x1> x1, T_sigma sigma, T_l l) {
  std::string message;
  try {
    stan::math::gp_exponential_cov(x1, sigma, l);
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

TEST(RevMath, gp_exp_quad_cov_vvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::gp_exponential_cov(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double l_v = l.val();
      double exp_val = exp(distance / (-l_v));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      // EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
      //     << "index: (" << i << ", " << j << ")";
      // EXPECT_FLOAT_EQ(
      //     sigma.val() * sigma.val() * exp_val * sq_distance / (sq_l * l.val()),
      //     grad[1])
      //     << "index: (" << i << ", " << j << ")";
      // EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * -distance / sq_l,
      //                 grad[2])
      //     << "index: (" << i << ", " << j << ")";
      // EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * distance / sq_l,
      //                 grad[3])
      //     << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}
