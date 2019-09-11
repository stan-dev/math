#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>
#include <string>
#include <vector>

template <typename T_x1, typename T_x2, typename T_sigma, typename T_l,
          typename T_p>
std::string pull_msg(std::vector<T_x1> x1, std::vector<T_x2> x2, T_sigma sigma,
                     T_l l, T_p p) {
  std::string message;
  try {
    stan::math::gp_periodic_cov(x1, x2, sigma, l, p);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

template <typename T_x1, typename T_sigma, typename T_l, typename T_p>
std::string pull_msg(std::vector<T_x1> x1, T_sigma sigma, T_l l, T_p p) {
  std::string message;
  try {
    stan::math::gp_periodic_cov(x1, sigma, l, p);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

TEST(RevMath, gp_periodic_cov_vvvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;
      stan::math::var p = 7;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      // Check positive definiteness
      Eigen::LLT<Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>>
          llt;
      llt.compute(cov);
      EXPECT_EQ(Eigen::ComputationInfo::Success, llt.info());

      // Check values
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val()
              * exp(-2.0
                    * pow(sin(M_PI * (x[i].val() - x[j].val()) / p.val()), 2)
                    / (l.val() * l.val())),
          cov(i, j).val())
          << "index: (" << i << ", " << j << ")";

      // Check gradients
      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(p);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);
      double distance = x[i].val() - x[j].val();
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_val_sq
                          / (sq_l * l.val()),
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI * distance / p.val() / p.val() / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * -4.0 * sin_cos_val
                          * M_PI / p.val() / sq_l,
                      grad[3])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI / p.val() / sq_l,
                      grad[4])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vvvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double p = 7;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p);
      double sin_cos_val = sin(M_PI * distance / p) * cos(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_val_sq
                          / (sq_l * l.val()),
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * -4.0 * sin_cos_val
                          * M_PI / p / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI / p / sq_l,
                      grad[3])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vvdv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      stan::math::var p = 7;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(p);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI * distance / p.val() / p.val() / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * -4.0 * sin_cos_val
                          * M_PI / p.val() / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI / p.val() / sq_l,
                      grad[3])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vdvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double sigma = 0.2;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var l = 5;
      stan::math::var p = 7;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);
      params.push_back(p);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_val_sq / (sq_l * l.val()),
          grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI
                          * distance / p.val() / p.val() / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * -4.0 * sin_cos_val * M_PI / p.val() / sq_l,
          grad[2])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI / p.val() / sq_l,
          grad[3])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vvdd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double l = 5;
  double p = 7;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p);
      double sin_cos_val = sin(M_PI * distance / p) * cos(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * -4.0 * sin_cos_val
                          * M_PI / p / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI / p / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vdvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double sigma = 0.2;
  double p = 7;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var l = 5;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p);
      double sin_cos_val = sin(M_PI * distance / p) * cos(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_val_sq / (sq_l * l.val()),
          grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * -4.0 * sin_cos_val * M_PI / p / sq_l,
          grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI / p / sq_l,
          grad[2])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vddv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double sigma = 0.2;
  double l = 5;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var p = 7;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(p);
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI
                          * distance / p.val() / p.val() / sq_l,
                      grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * -4.0 * sin_cos_val * M_PI / p.val() / sq_l,
          grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI / p.val() / sq_l,
          grad[2])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vddd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double sigma = 0.2;
  double l = 5;
  double p = 7;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(x[i]);
      params.push_back(x[j]);

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p);
      double sin_cos_val = sin(M_PI * distance / p) * cos(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * -4.0 * sin_cos_val * M_PI / p / sq_l,
          grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI / p / sq_l,
          grad[1])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_dvvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(p);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_val_sq
                          / (sq_l * l.val()),
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI * distance / p.val() / p.val() / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_dvvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double p = 7;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_val_sq
                          / (sq_l * l.val()),
                      grad[1])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}
//////////////////////////////////////////////////////
TEST(RevMath, gp_periodic_cov_dvdv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(p);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI * distance / p.val() / p.val() / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_ddvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double sigma = 0.2;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var l = 5;
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      // Check positive definiteness
      Eigen::LLT<Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>>
          llt;
      llt.compute(cov);
      EXPECT_EQ(Eigen::ComputationInfo::Success, llt.info());

      // Check values
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(-2.0 * pow(sin(M_PI * (x[i] - x[j]) / p.val()), 2)
                    / (l.val() * l.val())),
          cov(i, j).val())
          << "index: (" << i << ", " << j << ")";

      // Check gradients
      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);
      params.push_back(p);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_val_sq / (sq_l * l.val()),
          grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI
                          * distance / p.val() / p.val() / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_ddvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double sigma = 0.2;
  double p = 7;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_val_sq / (sq_l * l.val()),
          grad[0])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_dddv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double sigma = 0.2;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(p);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI
                          * distance / p.val() / p.val() / sq_l,
                      grad[0])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_dvdd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double l = 5;
  double p = 7;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);

      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_vvvv) {
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
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(p);
      params.push_back(x[i](0));
      params.push_back(x[i](1));
      params.push_back(x[j](0));
      params.push_back(x[j](1));

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_val_sq
                          / (sq_l * l.val()),
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI * distance / p.val() / p.val() / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";
      if (i == j) {
        EXPECT_FLOAT_EQ(0, grad[3]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[4]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[5]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[6]) << "index: (" << i << ", " << j << ")";
      } else {
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](0).val() - x[j](0).val()) / distance
                            / p.val() / sq_l,
                        grad[3])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](1).val() - x[j](1).val()) / distance
                            / p.val() / sq_l,
                        grad[4])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](0).val() - x[i](0).val()) / distance
                            / p.val() / sq_l,
                        grad[5])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](1).val() - x[i](1).val()) / distance
                            / p.val() / sq_l,
                        grad[6])
            << "index: (" << i << ", " << j << ")";
      }
      stan::math::recover_memory();
    }
  }
}

////////////////////////////////////////////////////////////////////////////
TEST(RevMath, gp_periodic_cov_vector_vvvd) {
  typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> vector_v;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double p = 7;
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

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

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
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p);
      double sin_cos_val = sin(M_PI * distance / p) * cos(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_val_sq
                          / (sq_l * l.val()),
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      if (i == j) {
        EXPECT_FLOAT_EQ(0, grad[2]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[3]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[4]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[5]) << "index: (" << i << ", " << j << ")";
      } else {
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](0).val() - x[j](0).val()) / distance / p
                            / sq_l,
                        grad[2])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](1).val() - x[j](1).val()) / distance / p
                            / sq_l,
                        grad[3])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](0).val() - x[i](0).val()) / distance / p
                            / sq_l,
                        grad[4])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](1).val() - x[i](1).val()) / distance / p
                            / sq_l,
                        grad[5])
            << "index: (" << i << ", " << j << ")";
      }
      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_vvdv) {
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
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(p);
      params.push_back(x[i](0));
      params.push_back(x[i](1));
      params.push_back(x[j](0));
      params.push_back(x[j](1));

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI * distance / p.val() / p.val() / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      if (i == j) {
        EXPECT_FLOAT_EQ(0, grad[2]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[3]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[4]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[5]) << "index: (" << i << ", " << j << ")";
      } else {
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](0).val() - x[j](0).val()) / distance
                            / p.val() / sq_l,
                        grad[2])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](1).val() - x[j](1).val()) / distance
                            / p.val() / sq_l,
                        grad[3])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](0).val() - x[i](0).val()) / distance
                            / p.val() / sq_l,
                        grad[4])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](1).val() - x[i](1).val()) / distance
                            / p.val() / sq_l,
                        grad[5])
            << "index: (" << i << ", " << j << ")";
      }
      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_vdvv) {
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
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);
      params.push_back(p);
      params.push_back(x[i](0));
      params.push_back(x[i](1));
      params.push_back(x[j](0));
      params.push_back(x[j](1));

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_val_sq / (sq_l * l.val()),
          grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI
                          * distance / p.val() / p.val() / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      if (i == j) {
        EXPECT_FLOAT_EQ(0, grad[2]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[3]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[4]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[5]) << "index: (" << i << ", " << j << ")";
      } else {
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](0).val() - x[j](0).val()) / distance
                            / p.val() / sq_l,
                        grad[2])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](1).val() - x[j](1).val()) / distance
                            / p.val() / sq_l,
                        grad[3])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](0).val() - x[i](0).val()) / distance
                            / p.val() / sq_l,
                        grad[4])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](1).val() - x[i](1).val()) / distance
                            / p.val() / sq_l,
                        grad[5])
            << "index: (" << i << ", " << j << ")";
      }
      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_vvdd) {
  typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> vector_v;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double l = 5;
  double p = 7;

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

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

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
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p);
      double sin_cos_val = sin(M_PI * distance / p) * cos(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      if (i == j) {
        EXPECT_FLOAT_EQ(0, grad[1]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[2]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[3]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[4]) << "index: (" << i << ", " << j << ")";
      } else {
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](0).val() - x[j](0).val()) / distance / p
                            / sq_l,
                        grad[1])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](1).val() - x[j](1).val()) / distance / p
                            / sq_l,
                        grad[2])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](0).val() - x[i](0).val()) / distance / p
                            / sq_l,
                        grad[3])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](1).val() - x[i](1).val()) / distance / p
                            / sq_l,
                        grad[4])
            << "index: (" << i << ", " << j << ")";
      }
      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_vdvd) {
  typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> vector_v;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double sigma = 0.2;
  double p = 7;

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

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

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
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p);
      double sin_cos_val = sin(M_PI * distance / p) * cos(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_val_sq / (sq_l * l.val()),
          grad[0])
          << "index: (" << i << ", " << j << ")";
      if (i == j) {
        EXPECT_FLOAT_EQ(0, grad[1]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[2]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[3]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[4]) << "index: (" << i << ", " << j << ")";
      } else {
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](0).val() - x[j](0).val()) / distance / p
                            / sq_l,
                        grad[1])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](1).val() - x[j](1).val()) / distance / p
                            / sq_l,
                        grad[2])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](0).val() - x[i](0).val()) / distance / p
                            / sq_l,
                        grad[3])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](1).val() - x[i](1).val()) / distance / p
                            / sq_l,
                        grad[4])
            << "index: (" << i << ", " << j << ")";
      }
      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_vddv) {
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

      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(p);
      params.push_back(x[i](0));
      params.push_back(x[i](1));
      params.push_back(x[j](0));
      params.push_back(x[j](1));

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI
                          * distance / p.val() / p.val() / sq_l,
                      grad[0])
          << "index: (" << i << ", " << j << ")";
      if (i == j) {
        EXPECT_FLOAT_EQ(0, grad[1]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[2]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[3]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[4]) << "index: (" << i << ", " << j << ")";
      } else {
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](0).val() - x[j](0).val()) / distance
                            / p.val() / sq_l,
                        grad[1])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](1).val() - x[j](1).val()) / distance
                            / p.val() / sq_l,
                        grad[2])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](0).val() - x[i](0).val()) / distance
                            / p.val() / sq_l,
                        grad[3])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](1).val() - x[i](1).val()) / distance
                            / p.val() / sq_l,
                        grad[4])
            << "index: (" << i << ", " << j << ")";
      }
      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_vddd) {
  typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> vector_v;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double sigma = 0.2;
  double l = 5;
  double p = 7;

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

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(x[i](0));
      params.push_back(x[i](1));
      params.push_back(x[j](0));
      params.push_back(x[j](1));

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p);
      double sin_cos_val = sin(M_PI * distance / p) * cos(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      if (i == j) {
        EXPECT_FLOAT_EQ(0, grad[0]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[1]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[2]) << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(0, grad[3]) << "index: (" << i << ", " << j << ")";
      } else {
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](0).val() - x[j](0).val()) / distance / p
                            / sq_l,
                        grad[0])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[i](1).val() - x[j](1).val()) / distance / p
                            / sq_l,
                        grad[1])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](0).val() - x[i](0).val()) / distance / p
                            / sq_l,
                        grad[2])
            << "index: (" << i << ", " << j << ")";
        EXPECT_FLOAT_EQ(-cov(i, j).val() * 4.0 * sin_cos_val * M_PI
                            * (x[j](1).val() - x[i](1).val()) / distance / p
                            / sq_l,
                        grad[3])
            << "index: (" << i << ", " << j << ")";
      }
      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_dvvv) {
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

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(p);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_val_sq
                          / (sq_l * l.val()),
                      grad[1])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI * distance / p.val() / p.val() / sq_l,
                      grad[2])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_dvvd) {
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

  double p = 7;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_val_sq
                          / (sq_l * l.val()),
                      grad[1])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_dvdv) {
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
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(p);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 4.0 * sin_cos_val
                          * M_PI * distance / p.val() / p.val() / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_ddvv) {
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
      stan::math::var p = 7;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);
      params.push_back(p);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_val_sq / (sq_l * l.val()),
          grad[0])
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI
                          * distance / p.val() / p.val() / sq_l,
                      grad[1])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_dvdd) {
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
  double p = 7;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);

      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val,
                      cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_ddvd) {
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
  double p = 5;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p);
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp_val * 4.0 * sin_val_sq / (sq_l * l.val()),
          grad[0])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov_vector_dddv) {
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
  double l = 5;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var p = 5;

      EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x, sigma, l, p));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(p);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_l = stan::math::square(l);
      double sin_val = sin(M_PI * distance / p.val());
      double sin_cos_val
          = sin(M_PI * distance / p.val()) * cos(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val * 4.0 * sin_cos_val * M_PI
                          * distance / p.val() / p.val() / sq_l,
                      grad[0])
          << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, gp_periodic_cov1_vec_eigen_rvec) {
  using stan::math::squared_distance;
  using stan::math::var;

  var sigma = 0.2;
  var l = 5;
  var p = 7;

  std::vector<Eigen::Matrix<var, 1, -1>> x1(3);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::Matrix<var, -1, -1> cov;
  EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x1, sigma, l, p));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(3, cov.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double distance = stan::math::distance(x1[i], x1[j]).val();
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);

      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(RevMath, gp_periodic_cov2_vec_eigen_rvec) {
  using stan::math::squared_distance;
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;
  var p = 7;

  std::vector<Eigen::Matrix<var, 1, -1>> x1(3);
  std::vector<Eigen::Matrix<var, 1, -1>> x2(4);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::Matrix<var, -1, -1> cov;
  EXPECT_NO_THROW(cov = stan::math::gp_periodic_cov(x1, x2, sigma, l, p));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      double distance = stan::math::distance(x1[i], x2[j]).val();
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);

      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val, cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::Matrix<var, -1, -1> cov2;
  cov2 = stan::math::gp_periodic_cov(x2, x1, sigma, l, p);
  EXPECT_EQ(4, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      double distance = stan::math::distance(x2[i], x1[j]).val();
      double sq_l = stan::math::square(l.val());
      double sin_val = sin(M_PI * distance / p.val());
      double sin_val_sq = stan::math::square(sin_val);
      double exp_val = exp(-2.0 * sin_val_sq / sq_l);

      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val, cov2(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov2(i, j).val(), cov(j, i).val());
    }
  }
}

TEST(RevMath, gp_periodic_cov2_vec_eigen_mixed) {
  using stan::math::squared_distance;
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;
  var p = 7;

  std::vector<Eigen::Matrix<var, 1, -1>> x1_rvec(3);
  std::vector<Eigen::Matrix<var, 1, -1>> x2_rvec(4);

  for (size_t i = 0; i < x1_rvec.size(); ++i) {
    x1_rvec[i].resize(1, 3);
    x1_rvec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_rvec.size(); ++i) {
    x2_rvec[i].resize(1, 3);
    x2_rvec[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<Eigen::Matrix<var, -1, 1>> x1_vec(3);
  std::vector<Eigen::Matrix<var, -1, 1>> x2_vec(4);

  for (size_t i = 0; i < x1_vec.size(); ++i) {
    x1_vec[i].resize(3, 1);
    x1_vec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_vec.size(); ++i) {
    x2_vec[i].resize(3, 1);
    x2_vec[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::Matrix<var, -1, -1> cov;
  EXPECT_NO_THROW(cov
                  = stan::math::gp_periodic_cov(x1_rvec, x2_vec, sigma, l, p));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val()
              * exp(-2.0
                    * stan::math::square(sin(
                          M_PI
                          * stan::math::distance(x1_rvec[i], x2_vec[j]).val()
                          / p.val()))
                    / l.val() / l.val()),
          cov(i, j).val())
          << "index: (" << i << ", " << j << ")";

  Eigen::Matrix<var, -1, -1> cov7;
  EXPECT_NO_THROW(cov7
                  = stan::math::gp_periodic_cov(x2_vec, x1_rvec, sigma, l, p));
  EXPECT_EQ(4, cov7.rows());
  EXPECT_EQ(3, cov7.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val()
              * exp(-2.0
                    * stan::math::square(sin(
                          M_PI
                          * stan::math::distance(x2_vec[i], x1_rvec[j]).val()
                          / p.val()))
                    / l.val() / l.val()),
          cov7(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov7(i, j).val(), cov(j, i).val());
    }

  Eigen::Matrix<var, -1, -1> cov2;
  EXPECT_NO_THROW(cov2
                  = stan::math::gp_periodic_cov(x1_vec, x2_rvec, sigma, l, p));
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(4, cov2.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val()
              * exp(-2.0
                    * stan::math::square(sin(
                          M_PI
                          * stan::math::distance(x1_vec[i], x2_rvec[j]).val()
                          / p.val()))
                    / l.val() / l.val()),
          cov2(i, j).val())
          << "index: (" << i << ", " << j << ")";

  Eigen::Matrix<var, -1, -1> cov8;
  EXPECT_NO_THROW(cov8
                  = stan::math::gp_periodic_cov(x2_rvec, x1_vec, sigma, l, p));
  EXPECT_EQ(4, cov8.rows());
  EXPECT_EQ(3, cov8.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val()
              * exp(-2.0
                    * stan::math::square(sin(
                          M_PI
                          * stan::math::distance(x2_rvec[i], x1_vec[j]).val()
                          / p.val()))
                    / l.val() / l.val()),
          cov8(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov8(i, j).val(), cov2(j, i).val());
    }

  Eigen::Matrix<var, -1, -1> cov3;
  EXPECT_NO_THROW(cov3
                  = stan::math::gp_periodic_cov(x2_vec, x2_rvec, sigma, l, p));
  EXPECT_EQ(4, cov3.rows());
  EXPECT_EQ(4, cov3.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val()
              * exp(-2.0
                    * stan::math::square(sin(
                          M_PI
                          * stan::math::distance(x2_vec[i], x2_rvec[j]).val()
                          / p.val()))
                    / l.val() / l.val()),
          cov3(i, j).val())
          << "index: (" << i << ", " << j << ")";

  Eigen::Matrix<var, -1, -1> cov4;
  EXPECT_NO_THROW(cov4
                  = stan::math::gp_periodic_cov(x2_rvec, x2_vec, sigma, l, p));
  EXPECT_EQ(4, cov4.rows());
  EXPECT_EQ(4, cov4.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val()
              * exp(-2.0
                    * stan::math::square(sin(
                          M_PI
                          * stan::math::distance(x2_rvec[i], x2_vec[j]).val()
                          / p.val()))
                    / l.val() / l.val()),
          cov4(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov4(i, j).val(), cov3(i, j).val());
    }

  Eigen::Matrix<var, -1, -1> cov5;
  EXPECT_NO_THROW(cov5
                  = stan::math::gp_periodic_cov(x1_rvec, x1_vec, sigma, l, p));
  EXPECT_EQ(3, cov5.rows());
  EXPECT_EQ(3, cov5.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val()
              * exp(-2.0
                    * stan::math::square(sin(
                          M_PI
                          * stan::math::distance(x1_rvec[i], x1_vec[j]).val()
                          / p.val()))
                    / l.val() / l.val()),
          cov5(i, j).val())
          << "index: (" << i << ", " << j << ")";

  Eigen::Matrix<var, -1, -1> cov6;
  EXPECT_NO_THROW(cov6
                  = stan::math::gp_periodic_cov(x1_vec, x1_rvec, sigma, l, p));
  EXPECT_EQ(3, cov6.rows());
  EXPECT_EQ(3, cov6.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma.val() * sigma.val()
              * exp(-2.0
                    * stan::math::square(sin(
                          M_PI
                          * stan::math::distance(x1_vec[i], x1_rvec[j]).val()
                          / p.val()))
                    / l.val() / l.val()),
          cov6(i, j).val())
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov6(i, j).val(), cov5(i, j).val());
    }
}

TEST(RevMath, gp_periodic_cov_domain_error_training) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;
  var p = 7;

  std::vector<var> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  std::vector<Eigen::Matrix<var, -1, 1>> x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1>> x_3(3);
  for (size_t i = 0; i < x_3.size(); ++i) {
    x_3[i].resize(1, 3);
    x_3[i] << 1, 2, 3;
  }

  var sigma_bad = -1;
  var l_bad = -1;
  var p_bad = -1;

  std::string msg1, msg2, msg3, msg4, msg5, msg6, msg7;
  msg1 = pull_msg(x, sigma, l, p_bad);
  msg2 = pull_msg(x, sigma, l_bad, p);
  msg3 = pull_msg(x, sigma_bad, l, p);
  msg4 = pull_msg(x, sigma, l_bad, p_bad);
  msg5 = pull_msg(x, sigma_bad, l, p_bad);
  msg6 = pull_msg(x, sigma_bad, l_bad, p);
  msg7 = pull_msg(x, sigma_bad, l_bad, p_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" period")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" length-scale")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" signal standard deviation"))
      << msg3;
  EXPECT_TRUE(std::string::npos != msg4.find(" length-scale")) << msg4;
  EXPECT_TRUE(std::string::npos != msg5.find(" signal standard deviation"))
      << msg5;
  EXPECT_TRUE(std::string::npos != msg6.find(" signal standard deviation"))
      << msg6;
  EXPECT_TRUE(std::string::npos != msg7.find(" signal standard deviation"))
      << msg7;

  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma_bad, l_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma_bad, l_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma_bad, l_bad, p_bad),
               std::domain_error);
}

TEST(RevMath, gp_periodic_cov_nan_error_training) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;
  var p = 5;

  std::vector<var> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  std::vector<Eigen::Matrix<var, -1, 1>> x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1>> x_3(3);
  for (size_t i = 0; i < x_3.size(); ++i) {
    x_3[i].resize(1, 3);
    x_3[i] << 1, 2, 3;
  }

  std::vector<var> x_bad(x);
  x_bad[1] = std::numeric_limits<var>::quiet_NaN();

  std::vector<Eigen::Matrix<var, -1, 1>> x_bad_2(x_2);
  x_bad_2[1](1) = std::numeric_limits<var>::quiet_NaN();

  std::vector<Eigen::Matrix<var, 1, -1>> x_bad_3(x_3);
  x_bad_3[1](1) = std::numeric_limits<var>::quiet_NaN();

  var sigma_bad = std::numeric_limits<var>::quiet_NaN();
  var l_bad = std::numeric_limits<var>::quiet_NaN();
  var p_bad = std::numeric_limits<var>::quiet_NaN();

  std::string msg1, msg2, msg3, msg4, msg5, msg6, msg7;
  msg1 = pull_msg(x, sigma, l, p_bad);
  msg2 = pull_msg(x, sigma, l_bad, p);
  msg3 = pull_msg(x, sigma_bad, l, p);
  msg4 = pull_msg(x, sigma, l_bad, p_bad);
  msg5 = pull_msg(x, sigma_bad, l, p_bad);
  msg6 = pull_msg(x, sigma_bad, l_bad, p);
  msg7 = pull_msg(x, sigma_bad, l_bad, p_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" period")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" length-scale")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" signal standard deviation"))
      << msg3;
  EXPECT_TRUE(std::string::npos != msg4.find(" length-scale")) << msg4;
  EXPECT_TRUE(std::string::npos != msg5.find(" signal standard deviation"))
      << msg5;
  EXPECT_TRUE(std::string::npos != msg6.find(" signal standard deviation"))
      << msg6;
  EXPECT_TRUE(std::string::npos != msg7.find(" signal standard deviation"))
      << msg7;

  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_2, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_2, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_2, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_2, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_2, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_2, sigma_bad, l_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_3, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_3, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_3, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_3, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_3, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_3, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_3, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_3, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_bad_3, sigma_bad, l_bad, p_bad),
               std::domain_error);
}

TEST(RevMath, gp_periodic_cov_domain_error) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;
  var p = 7;

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
  var p_bad = -1;

  std::string msg1, msg2, msg3, msg4, msg5, msg6, msg7;
  msg1 = pull_msg(x1, x2, sigma, l, p_bad);
  msg2 = pull_msg(x1, x2, sigma, l_bad, p);
  msg3 = pull_msg(x1, x2, sigma_bad, l, p);
  msg4 = pull_msg(x1, x2, sigma, l_bad, p_bad);
  msg5 = pull_msg(x1, x2, sigma_bad, l, p_bad);
  msg6 = pull_msg(x1, x2, sigma_bad, l_bad, p);
  msg7 = pull_msg(x1, x2, sigma_bad, l_bad, p_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" period")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" length-scale")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" signal standard deviation"))
      << msg3;
  EXPECT_TRUE(std::string::npos != msg4.find(" length-scale")) << msg4;
  EXPECT_TRUE(std::string::npos != msg5.find(" signal standard deviation"))
      << msg5;
  EXPECT_TRUE(std::string::npos != msg6.find(" signal standard deviation"))
      << msg6;
  EXPECT_TRUE(std::string::npos != msg7.find(" signal standard deviation"))
      << msg7;

  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma_bad, l_bad, p_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma_bad, l_bad, p_bad),
      std::domain_error);

  std::vector<Eigen::Matrix<var, 1, -1>> x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1>> x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 3);
    x_rvec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma_bad, l_bad, p_bad),
      std::domain_error);
}

TEST(RevMath, gp_periodic_cov2_nan_domain_error) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;
  var p = 7;

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
  var p_bad = std::numeric_limits<var>::quiet_NaN();

  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2, sigma_bad, l_bad, p_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma_bad, l_bad, p_bad),
      std::domain_error);

  std::vector<Eigen::Matrix<var, 1, -1>> x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1>> x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 3);
    x_rvec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1, x_rvec_2, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1, x_vec_2, sigma_bad, l_bad, p_bad),
      std::domain_error);

  std::vector<var> x1_bad(x1);
  x1_bad[1] = std::numeric_limits<var>::quiet_NaN();
  std::vector<var> x2_bad(x2);
  x2_bad[1] = std::numeric_limits<var>::quiet_NaN();

  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_1_bad(x_vec_1);
  x_vec_1_bad[1](1) = std::numeric_limits<var>::quiet_NaN();
  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_2_bad(x_vec_2);
  x_vec_2_bad[1](1) = std::numeric_limits<var>::quiet_NaN();

  std::vector<Eigen::Matrix<var, 1, -1>> x_rvec_1_bad(x_rvec_1);
  x_rvec_1_bad[1](1) = std::numeric_limits<var>::quiet_NaN();
  std::vector<Eigen::Matrix<var, 1, -1>> x_rvec_2_bad(x_rvec_2);
  x_rvec_2_bad[1](1) = std::numeric_limits<var>::quiet_NaN();

  EXPECT_THROW(stan::math::gp_periodic_cov(x1_bad, x2, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1, x2_bad, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x1_bad, x2_bad, sigma, l, p),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1_bad, x_vec_2, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_vec_2_bad, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1_bad, x_vec_2_bad, sigma, l, p),
      std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1_bad, x_rvec_2, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_rvec_2_bad, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1_bad, x_rvec_2_bad, sigma, l, p),
      std::domain_error);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1_bad, x_rvec_2, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1_bad, x_vec_2, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_rvec_2_bad, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_vec_2_bad, sigma, l, p),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_vec_1_bad, x_rvec_2_bad, sigma, l, p),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gp_periodic_cov(x_rvec_1_bad, x_vec_2_bad, sigma, l, p),
      std::domain_error);
}

TEST(RevMath, gp_periodic_cov2_dim_mismatch_vec_eigen_vec) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;
  var p = 7;

  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(4, 1);
    x_vec_2[i] << 4, 1, 3, 1;
  }

  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma, l, p),
               std::invalid_argument);
}

TEST(RevMath, gp_periodic_cov2_dim_mismatch_vec_eigen_rvec) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;
  var p = 7;

  std::vector<Eigen::Matrix<var, 1, -1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(1, 3);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(1, 4);
    x_vec_2[i] << 4, 1, 3, 1;
  }

  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_vec_2, sigma, l, p),
               std::invalid_argument);
}

TEST(RevMath, gp_periodic_cov2_dim_mismatch_vec_eigen_mixed) {
  using stan::math::var;
  var sigma = 0.2;
  var l = 5;
  var p = 7;

  std::vector<Eigen::Matrix<var, 1, -1>> x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_1(4);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(4, 1);
    x_vec_1[i] << 4, 1, 3, 1;
  }

  std::vector<Eigen::Matrix<var, -1, 1>> x_vec_2(3);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<var, 1, -1>> x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 4);
    x_rvec_2[i] << 1, 2, 3, 4;
  }

  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_1, x_vec_1, sigma, l, p),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_1, x_rvec_1, sigma, l, p),
               std::invalid_argument);

  EXPECT_THROW(stan::math::gp_periodic_cov(x_rvec_2, x_vec_2, sigma, l, p),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gp_periodic_cov(x_vec_2, x_rvec_2, sigma, l, p),
               std::invalid_argument);
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::to_var;
  std::vector<double> x(3);
  double sigma = 0.2;
  double l = 5;
  double p = 5;
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  test::check_varis_on_stack(stan::math::gp_periodic_cov(
      to_var(x), to_var(sigma), to_var(l), to_var(p)));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(to_var(x), to_var(sigma), to_var(l), p));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(to_var(x), to_var(sigma), l, to_var(p)));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(to_var(x), sigma, to_var(l), to_var(p)));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(to_var(x), to_var(sigma), l, p));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(to_var(x), sigma, to_var(l), p));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(to_var(x), sigma, l, to_var(p)));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(to_var(x), sigma, l, p));

  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(x, to_var(sigma), to_var(l), to_var(p)));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(x, to_var(sigma), to_var(l), p));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(x, to_var(sigma), l, to_var(p)));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(x, sigma, to_var(l), to_var(p)));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(x, to_var(sigma), l, p));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(x, sigma, to_var(l), p));
  test::check_varis_on_stack(
      stan::math::gp_periodic_cov(x, sigma, l, to_var(p)));
}
