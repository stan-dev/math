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
    stan::math::cov_exp_quad(x1, x2, sigma, l);
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
    stan::math::cov_exp_quad(x1, sigma, l);
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

TEST(RevMath, cov_exp_quad_vvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val() * exp_val * sq_distance / (sq_l * l.val()),
          grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * -distance / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * distance / sq_l,
                      grad[3])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * -distance / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * distance / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vdv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double sigma = 0.2;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var l = 5;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * sq_distance / (sq_l * l.val()),
                      grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * -distance / sq_l, grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * distance / sq_l, grad[2])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vdd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double sigma = 0.2;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * -distance / sq_l, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * distance / sq_l, grad[1])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_dvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val() * exp_val * sq_distance / (sq_l * l.val()),
          grad[1])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_dvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_ddv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double sigma = 0.2;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * sq_distance / (sq_l * l.val()),
                      grad[0]);

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vector_vvv) {
  typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> vector_v;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<vector_v> x(3);
      vector_v x0(2), x1(2), x2(2);
      x0(0) = -2;
      x0(1) = -2;
      x1(0) = 1;
      x1(1) = 2;
      x2(0) = -0.5;
      x2(1) = 0.0;

      x[0] = x0;
      x[1] = x1;
      x[2] = x2;

      stan::math::var sigma = 0.2;
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(x[i](0));
      params.push_back(x[i](1));
      params.push_back(x[j](0));
      params.push_back(x[j](1));

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val() * exp_val * sq_distance / (sq_l * l.val()),
          grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[i](0).val() - x[j](0).val()) / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[i](1).val() - x[j](1).val()) / sq_l,
                      grad[3])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[j](0).val() - x[i](0).val()) / sq_l,
                      grad[4])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[j](1).val() - x[i](1).val()) / sq_l,
                      grad[5])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vector_vvd) {
  typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> vector_v;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<vector_v> x(3);
      vector_v x0(2), x1(2), x2(2);
      x0(0) = -2;
      x0(1) = -2;
      x1(0) = 1;
      x1(1) = 2;
      x2(0) = -0.5;
      x2(1) = 0.0;

      x[0] = x0;
      x[1] = x1;
      x[2] = x2;

      stan::math::var sigma = 0.2;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(x[i](0));
      params.push_back(x[i](1));
      params.push_back(x[j](0));
      params.push_back(x[j](1));

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[i](0).val() - x[j](0).val()) / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[i](1).val() - x[j](1).val()) / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[j](0).val() - x[i](0).val()) / sq_l,
                      grad[3])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[j](1).val() - x[i](1).val()) / sq_l,
                      grad[4])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vector_vdv) {
  typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> vector_v;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double sigma = 0.2;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<vector_v> x(3);
      vector_v x0(2), x1(2), x2(2);
      x0(0) = -2;
      x0(1) = -2;
      x1(0) = 1;
      x1(1) = 2;
      x2(0) = -0.5;
      x2(1) = 0.0;

      x[0] = x0;
      x[1] = x1;
      x[2] = x2;

      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);
      params.push_back(x[i](0));
      params.push_back(x[i](1));
      params.push_back(x[j](0));
      params.push_back(x[j](1));

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * sq_distance / (sq_l * l.val()),
                      grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[i](0).val() - x[j](0).val()) / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[i](1).val() - x[j](1).val()) / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[j](0).val() - x[i](0).val()) / sq_l,
                      grad[3])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[j](1).val() - x[i](1).val()) / sq_l,
                      grad[4])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vector_vdd) {
  typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> vector_v;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double sigma = 0.2;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<vector_v> x(3);
      vector_v x0(2), x1(2), x2(2);
      x0(0) = -2;
      x0(1) = -2;
      x1(0) = 1;
      x1(1) = 2;
      x2(0) = -0.5;
      x2(1) = 0.0;

      x[0] = x0;
      x[1] = x1;
      x[2] = x2;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(x[i](0));
      params.push_back(x[i](1));
      params.push_back(x[j](0));
      params.push_back(x[j](1));

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[i](0).val() - x[j](0).val()) / sq_l,
                      grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[i](1).val() - x[j](1).val()) / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[j](0).val() - x[i](0).val()) / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val() * (x[j](1).val() - x[i](1).val()) / sq_l,
                      grad[3])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vector_dvv) {
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_d;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  std::vector<vector_d> x(3);
  vector_d x0(2), x1(2), x2(2);
  x0(0) = -2;
  x0(1) = -2;
  x1(0) = 1;
  x1(1) = 2;
  x2(0) = -0.5;
  x2(1) = 0.0;

  x[0] = x0;
  x[1] = x1;
  x[2] = x2;

  // std::cout << "x[0]: " << x[0] << std::endl
  //           << "x[1]: " << x[1] << std::endl
  //           << "x[2]: " << x[2] << std::endl
  //           << std::endl;
  // std::cout << "x[0]: " << x[0] << std::endl
  //           << "x[1]: " << x[1] << std::endl
  //           << "x[2]: " << x[2] << std::endl
  //           << std::endl;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);

      cov(i, j).grad(params, grad);

      // std::cout << "x[0]: " << x[0] << std::endl
      //           << "x[1]: " << x[1] << std::endl
      //           << "x[2]: " << x[2] << std::endl
      //           << std::endl;
      double sq_distance = stan::math::squared_distance(
          stan::math::value_of(x[i]), stan::math::value_of(x[j]));

      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      // << std::endl
      // << "x[i]: " << x[i] << std::endl
      // << "x[j]: " << x[j] << std::endl
      // << "dist: " << distance << std::endl
      // << "sq_dist: " << sq_distance << std::endl
      // << "sq_l: " << sq_l << std::endl
      // << "exp_val: " << exp_val << std::endl;
      // EXPECT_FLOAT_EQ(- 2 * cov(i, j).val() / pow(l.val(), 3),
      //                 grad[1])
      //   << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vector_dvd) {
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_d;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<vector_d> x(3);
  vector_d x0(2), x1(2), x2(2);
  x0(0) = -2;
  x0(1) = -2;
  x1(0) = 1;
  x1(1) = 2;
  x2(0) = -0.5;
  x2(1) = 0.0;

  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad_vector_ddv) {
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_d;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<vector_d> x(3);
  vector_d x0(2), x1(2), x2(2);
  x0(0) = -2;
  x0(1) = -2;
  x1(0) = 1;
  x1(1) = 2;
  x2(0) = -0.5;
  x2(1) = 0.0;

  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
  double sigma = 0.2;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * sq_distance / (sq_l * l.val()),
                      grad[0])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_exp_quad1_vec_eigen_rvec) {
  using stan::math::squared_distance;
  using stan::math::var;

  var sigma = 0.2;
  var l = 5;

  std::vector<Eigen::Matrix<var, 1, -1> > x1(3);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::Matrix<var, -1, -1> cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(3, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x1[i], x1[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
}

TEST(RevMath, cov_exp_quad2_vec_eigen_rvec) {
  using stan::math::squared_distance;
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;

  std::vector<Eigen::Matrix<var, 1, -1> > x1(3);
  std::vector<Eigen::Matrix<var, 1, -1> > x2(4);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::Matrix<var, -1, -1> cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x1[i], x2[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";

  Eigen::Matrix<var, -1, -1> cov2;
  cov2 = stan::math::cov_exp_quad(x2, x1, sigma, l);
  EXPECT_EQ(4, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x2[i], x1[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov2(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov2(i, j).val(), cov(j, i).val());
    }
}

TEST(RevMath, cov_exp_quad2_vec_eigen_mixed) {
  using stan::math::squared_distance;
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;

  std::vector<Eigen::Matrix<var, 1, -1> > x1_rvec(3);
  std::vector<Eigen::Matrix<var, 1, -1> > x2_rvec(4);

  for (size_t i = 0; i < x1_rvec.size(); ++i) {
    x1_rvec[i].resize(1, 3);
    x1_rvec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_rvec.size(); ++i) {
    x2_rvec[i].resize(1, 3);
    x2_rvec[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<Eigen::Matrix<var, -1, 1> > x1_vec(3);
  std::vector<Eigen::Matrix<var, -1, 1> > x2_vec(4);

  for (size_t i = 0; i < x1_vec.size(); ++i) {
    x1_vec[i].resize(3, 1);
    x1_vec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_vec.size(); ++i) {
    x2_vec[i].resize(3, 1);
    x2_vec[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::Matrix<var, -1, -1> cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1_rvec, x2_vec, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x1_rvec[i], x2_vec[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";

  Eigen::Matrix<var, -1, -1> cov7;
  EXPECT_NO_THROW(cov7 = stan::math::cov_exp_quad(x2_vec, x1_rvec, sigma, l));
  EXPECT_EQ(4, cov7.rows());
  EXPECT_EQ(3, cov7.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x2_vec[i], x1_rvec[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov7(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov7(i, j).val(), cov(j, i).val());
    }

  Eigen::Matrix<var, -1, -1> cov2;
  EXPECT_NO_THROW(cov2 = stan::math::cov_exp_quad(x1_vec, x2_rvec, sigma, l));
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(4, cov2.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x1_vec[i], x2_rvec[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov2(i, j).val())
          << "index: (" << i << ", " << j << ")";

  Eigen::Matrix<var, -1, -1> cov8;
  EXPECT_NO_THROW(cov8 = stan::math::cov_exp_quad(x2_rvec, x1_vec, sigma, l));
  EXPECT_EQ(4, cov8.rows());
  EXPECT_EQ(3, cov8.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x2_rvec[i], x1_vec[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov8(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov8(i, j).val(), cov2(j, i).val());
    }

  Eigen::Matrix<var, -1, -1> cov3;
  EXPECT_NO_THROW(cov3 = stan::math::cov_exp_quad(x2_vec, x2_rvec, sigma, l));
  EXPECT_EQ(4, cov3.rows());
  EXPECT_EQ(4, cov3.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x2_vec[i], x2_rvec[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov3(i, j).val())
          << "index: (" << i << ", " << j << ")";

  Eigen::Matrix<var, -1, -1> cov4;
  EXPECT_NO_THROW(cov4 = stan::math::cov_exp_quad(x2_rvec, x2_vec, sigma, l));
  EXPECT_EQ(4, cov4.rows());
  EXPECT_EQ(4, cov4.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x2_rvec[i], x2_vec[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov4(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov4(i, j).val(), cov3(i, j).val());
    }

  Eigen::Matrix<var, -1, -1> cov5;
  EXPECT_NO_THROW(cov5 = stan::math::cov_exp_quad(x1_rvec, x1_vec, sigma, l));
  EXPECT_EQ(3, cov5.rows());
  EXPECT_EQ(3, cov5.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x1_rvec[i], x1_vec[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov5(i, j).val())
          << "index: (" << i << ", " << j << ")";

  Eigen::Matrix<var, -1, -1> cov6;
  EXPECT_NO_THROW(cov6 = stan::math::cov_exp_quad(x1_vec, x1_rvec, sigma, l));
  EXPECT_EQ(3, cov6.rows());
  EXPECT_EQ(3, cov6.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val()
                          * exp(squared_distance(x1_vec[i], x1_rvec[j]).val()
                                / (-2.0 * l.val() * l.val())),
                      cov6(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov6(i, j).val(), cov5(i, j).val());
    }
}

TEST(RevMath, cov_exp_quad_domain_error_training) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;

  std::vector<var> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  std::vector<Eigen::Matrix<var, -1, 1> > x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1> > x_3(3);
  for (size_t i = 0; i < x_3.size(); ++i) {
    x_3[i].resize(1, 3);
    x_3[i] << 1, 2, 3;
  }

  var sigma_bad = -1;
  var l_bad = -1;

  std::string msg1, msg2, msg3;
  msg1 = pull_msg(x, sigma, l_bad);
  msg2 = pull_msg(x, sigma_bad, l);
  msg3 = pull_msg(x, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" magnitude")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" magnitude")) << msg3;

  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l_bad),
               std::domain_error);
}

TEST(RevMath, cov_exp_quad_nan_error_training) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;

  std::vector<var> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  std::vector<Eigen::Matrix<var, -1, 1> > x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1> > x_3(3);
  for (size_t i = 0; i < x_3.size(); ++i) {
    x_3[i].resize(1, 3);
    x_3[i] << 1, 2, 3;
  }

  std::vector<var> x_bad(x);
  x_bad[1] = std::numeric_limits<var>::quiet_NaN();

  std::vector<Eigen::Matrix<var, -1, 1> > x_bad_2(x_2);
  x_bad_2[1](1) = std::numeric_limits<var>::quiet_NaN();

  std::vector<Eigen::Matrix<var, 1, -1> > x_bad_3(x_3);
  x_bad_3[1](1) = std::numeric_limits<var>::quiet_NaN();

  var sigma_bad = std::numeric_limits<var>::quiet_NaN();
  var l_bad = std::numeric_limits<var>::quiet_NaN();

  std::string msg1, msg2, msg3;
  msg1 = pull_msg(x, sigma, l_bad);
  msg2 = pull_msg(x, sigma_bad, l);
  msg3 = pull_msg(x, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" magnitude")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" magnitude")) << msg3;

  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad, sigma, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma_bad, l_bad),
               std::domain_error);
}

TEST(RevMath, cov_exp_quad_domain_error) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;

  std::vector<var> x1(3);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  std::vector<var> x2(4);
  x2[0] = -2;
  x2[1] = -1;
  x2[2] = -0.5;
  x2[3] = -5;

  var sigma_bad = -1;
  var l_bad = -1;

  std::string msg1, msg2, msg3;
  msg1 = pull_msg(x1, x2, sigma, l_bad);
  msg2 = pull_msg(x1, x2, sigma_bad, l);
  msg3 = pull_msg(x1, x2, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" magnitude")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" magnitude")) << msg3;

  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<var, 1, -1> > x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1> > x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 3);
    x_rvec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l_bad),
               std::domain_error);
}

TEST(RevMath, cov_exp_quad2_nan_domain_error) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;

  std::vector<var> x1(3);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  std::vector<var> x2(4);
  x2[0] = -2;
  x2[1] = -1;
  x2[2] = -0.5;
  x2[3] = -5;

  var sigma_bad = std::numeric_limits<var>::quiet_NaN();
  var l_bad = std::numeric_limits<var>::quiet_NaN();

  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<var, 1, -1> > x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1> > x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 3);
    x_rvec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<var> x1_bad(x1);
  x1_bad[1] = std::numeric_limits<var>::quiet_NaN();
  std::vector<var> x2_bad(x2);
  x2_bad[1] = std::numeric_limits<var>::quiet_NaN();

  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_1_bad(x_vec_1);
  x_vec_1_bad[1](1) = std::numeric_limits<var>::quiet_NaN();
  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_2_bad(x_vec_2);
  x_vec_2_bad[1](1) = std::numeric_limits<var>::quiet_NaN();

  std::vector<Eigen::Matrix<var, 1, -1> > x_rvec_1_bad(x_rvec_1);
  x_rvec_1_bad[1](1) = std::numeric_limits<var>::quiet_NaN();
  std::vector<Eigen::Matrix<var, 1, -1> > x_rvec_2_bad(x_rvec_2);
  x_rvec_2_bad[1](1) = std::numeric_limits<var>::quiet_NaN();

  EXPECT_THROW(stan::math::cov_exp_quad(x1_bad, x2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1_bad, x2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_vec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_vec_2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_rvec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_rvec_2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_rvec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_vec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_vec_2_bad, sigma, l),
               std::domain_error);
}

TEST(RevMath, cov_exp_quad2_dim_mismatch_vec_eigen_vec) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;

  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(4, 1);
    x_vec_2[i] << 4, 1, 3, 1;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(RevMath, cov_exp_quad2_dim_mismatch_vec_eigen_rvec) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;

  std::vector<Eigen::Matrix<var, 1, -1> > x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(1, 3);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1> > x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(1, 4);
    x_vec_2[i] << 4, 1, 3, 1;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(RevMath, cov_exp_quad2_dim_mismatch_vec_eigen_mixed) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;

  std::vector<Eigen::Matrix<var, 1, -1> > x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_1(4);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(4, 1);
    x_vec_1[i] << 4, 1, 3, 1;
  }

  std::vector<Eigen::Matrix<var, -1, 1> > x_vec_2(3);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1> > x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 4);
    x_rvec_2[i] << 1, 2, 3, 4;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_1, sigma, l),
               std::invalid_argument);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_1, sigma, l),
               std::invalid_argument);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_2, x_vec_2, sigma, l),
               std::invalid_argument);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_2, x_rvec_2, sigma, l),
               std::invalid_argument);
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::to_var;
  std::vector<double> x(3);
  double sigma = 0.2;
  double l = 5;
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  test::check_varis_on_stack(
      stan::math::cov_exp_quad(to_var(x), to_var(sigma), to_var(l)));
  test::check_varis_on_stack(
      stan::math::cov_exp_quad(to_var(x), to_var(sigma), l));
  test::check_varis_on_stack(
      stan::math::cov_exp_quad(to_var(x), sigma, to_var(l)));
  test::check_varis_on_stack(stan::math::cov_exp_quad(to_var(x), sigma, l));
  test::check_varis_on_stack(
      stan::math::cov_exp_quad(x, to_var(sigma), to_var(l)));
  test::check_varis_on_stack(stan::math::cov_exp_quad(x, to_var(sigma), l));
  test::check_varis_on_stack(stan::math::cov_exp_quad(x, sigma, to_var(l)));
}
