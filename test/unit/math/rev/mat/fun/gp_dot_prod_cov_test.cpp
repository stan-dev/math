#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

template <typename T_x, typename T_sigma>
std::string pull_msg(std::vector<T_x> x, T_sigma sigma) {
  std::string message;
  try {
    stan::math::gp_dot_prod_cov(x, sigma);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exection";
  }
  return message;
}

template <typename T_x1, typename T_x2, typename T_sigma>
std::string pull_msg(std::vector<T_x1> x1, std::vector<T_x2> x2,
                     T_sigma sigma) {
  std::string message;
  try {
    stan::math::gp_dot_prod_cov(x1, x2, sigma);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exection";
  }
  return message;
}

TEST(RevMath, gp_dot_prod_cov_vvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;

      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;      

      //EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma));
      
      // std::vector<double> grad;
      // std::vector<stan::math::var> params;
      // params.push_back(sigma);
      // params.push_back(x[i]);
      // params.push_back(x[j]);


      // cov(i, j).grad(params, grad);

    }
  }
}

TEST(RevMath, gpdot_prod_cov_nan_error_training) {
  using stan::math::var;
  var sigma = 0.2;

  std::vector<var> x_0(3);
  x_0[0] = -2;
  x_0[1] = -1;
  x_0[2] = -0.5;

  std::vector<Eigen::Matrix<var, -1, 1> > x_1(3);
  for (size_t i = 0; i < x_1.size(); ++i) {
    x_1[i].resize(3, 1);
    x_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1> > x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(1, 3);
    x_2[i] << 1, 2, 3;
  }

  Eigen::Matrix<var, -1, -1> cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x_1, x_2, sigma));  // function call works
  
}
