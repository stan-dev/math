#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

TEST(MathPrimMat, vec_double_cov_exp_quad1) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd cov;
  cov = stan::math::cov_exp_quad(x, sigma, l);
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma * sigma * exp(pow(x[i] - x[j], 2) / (-2.0 * l * l)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_cov_exp_quad1) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, -1, 1> > x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::cov_exp_quad(x, sigma, l);
  EXPECT_NO_THROW(stan::math::cov_exp_quad(x, sigma, l));
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x[i], x[j]) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, rvec_eigen_cov_exp_quad1) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1> > x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(1, 3);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::cov_exp_quad(x, sigma, l);
  EXPECT_NO_THROW(stan::math::cov_exp_quad(x, sigma, l));
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x[i], x[j]) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_double_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x1(3);
  std::vector<double> x2(4);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  x2[0] = 5;
  x2[1] = 0;
  x2[2] = -4;
  x2[3] = 1.1;

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(pow(x1[i] - x2[j], 2) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_rvec_cov_exp_quad2) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1> > x1(3);
  std::vector<Eigen::Matrix<double, 1, -1> > x2(4);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x1[i], x2[j]) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov2;
  cov2 = stan::math::cov_exp_quad(x2, x1, sigma, l);
  EXPECT_EQ(4, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x2[i], x1[j]) / (-2.0 * l * l)),
          cov2(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov2(i, j), cov(j, i));
    }
}

TEST(MathPrimMat, vec_eigen_vec_cov_exp_quad2) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, -1, 1> > x1(3);
  std::vector<Eigen::Matrix<double, -1, 1> > x2(4);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(3, 1);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x1[i], x2[j]) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov2;
  cov2 = stan::math::cov_exp_quad(x2, x1, sigma, l);
  EXPECT_EQ(4, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x2[i], x1[j]) / (-2.0 * l * l)),
          cov2(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov2(i, j), cov(j, i));
    }
}

TEST(MathPrimMat, vec_eigen_mixed_cov_exp_quad2) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1> > x1_rvec(3);
  std::vector<Eigen::Matrix<double, 1, -1> > x2_rvec(4);

  for (size_t i = 0; i < x1_rvec.size(); ++i) {
    x1_rvec[i].resize(1, 3);
    x1_rvec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_rvec.size(); ++i) {
    x2_rvec[i].resize(1, 3);
    x2_rvec[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x1_vec(3);
  std::vector<Eigen::Matrix<double, -1, 1> > x2_vec(4);

  for (size_t i = 0; i < x1_vec.size(); ++i) {
    x1_vec[i].resize(3, 1);
    x1_vec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_vec.size(); ++i) {
    x2_vec[i].resize(3, 1);
    x2_vec[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1_rvec, x2_vec, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(squared_distance(x1_rvec[i], x2_vec[j]) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov7;
  EXPECT_NO_THROW(cov7 = stan::math::cov_exp_quad(x2_vec, x1_rvec, sigma, l));
  EXPECT_EQ(4, cov7.rows());
  EXPECT_EQ(3, cov7.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(squared_distance(x2_vec[i], x1_rvec[j]) / (-2.0 * l * l)),
          cov7(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov7(i, j), cov(j, i));
    }

  Eigen::MatrixXd cov2;
  EXPECT_NO_THROW(cov2 = stan::math::cov_exp_quad(x1_vec, x2_rvec, sigma, l));
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(4, cov2.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(squared_distance(x1_vec[i], x2_rvec[j]) / (-2.0 * l * l)),
          cov2(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov8;
  EXPECT_NO_THROW(cov8 = stan::math::cov_exp_quad(x2_rvec, x1_vec, sigma, l));
  EXPECT_EQ(4, cov8.rows());
  EXPECT_EQ(3, cov8.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(squared_distance(x2_rvec[i], x1_vec[j]) / (-2.0 * l * l)),
          cov8(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov8(i, j), cov2(j, i));
    }

  Eigen::MatrixXd cov3;
  EXPECT_NO_THROW(cov3 = stan::math::cov_exp_quad(x2_vec, x2_rvec, sigma, l));
  EXPECT_EQ(4, cov3.rows());
  EXPECT_EQ(4, cov3.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(squared_distance(x2_vec[i], x2_rvec[j]) / (-2.0 * l * l)),
          cov3(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov4;
  EXPECT_NO_THROW(cov4 = stan::math::cov_exp_quad(x2_rvec, x2_vec, sigma, l));
  EXPECT_EQ(4, cov4.rows());
  EXPECT_EQ(4, cov4.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(squared_distance(x2_rvec[i], x2_vec[j]) / (-2.0 * l * l)),
          cov4(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov4(i, j), cov3(i, j));
    }

  Eigen::MatrixXd cov5;
  EXPECT_NO_THROW(cov5 = stan::math::cov_exp_quad(x1_rvec, x1_vec, sigma, l));
  EXPECT_EQ(3, cov5.rows());
  EXPECT_EQ(3, cov5.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(squared_distance(x1_rvec[i], x1_vec[j]) / (-2.0 * l * l)),
          cov5(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov6;
  EXPECT_NO_THROW(cov6 = stan::math::cov_exp_quad(x1_vec, x1_rvec, sigma, l));
  EXPECT_EQ(3, cov6.rows());
  EXPECT_EQ(3, cov6.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(squared_distance(x1_vec[i], x1_rvec[j]) / (-2.0 * l * l)),
          cov6(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov6(i, j), cov5(i, j));
    }
}

TEST(MathPrimMat, domain_error_training_sig_l_cov_exp_quad) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  std::vector<Eigen::Matrix<double, -1, 1> > x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x_3(3);
  for (size_t i = 0; i < x_3.size(); ++i) {
    x_3[i].resize(1, 3);
    x_3[i] << 1, 2, 3;
  }

  double sigma_bad = -1;
  double l_bad = -1;

  EXPECT_THROW_MSG(stan::math::cov_exp_quad(x, sigma, l_bad), std::domain_error,
                   " length scale");
  EXPECT_THROW_MSG(stan::math::cov_exp_quad(x, sigma_bad, l), std::domain_error,
                   " magnitude");
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

TEST(MathPrimMat, nan_error_training_sig_l_cov_exp_quad) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  std::vector<Eigen::Matrix<double, -1, 1> > x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x_3(3);
  for (size_t i = 0; i < x_3.size(); ++i) {
    x_3[i].resize(1, 3);
    x_3[i] << 1, 2, 3;
  }

  std::vector<double> x_bad(x);
  x_bad[1] = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1> > x_bad_2(x_2);
  x_bad_2[1](1) = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, 1, -1> > x_bad_3(x_3);
  x_bad_3[1](1) = std::numeric_limits<double>::quiet_NaN();

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();
  double l_bad = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW_MSG(stan::math::cov_exp_quad(x, sigma, l_bad), std::domain_error,
                   " length scale");
  EXPECT_THROW_MSG(stan::math::cov_exp_quad(x, sigma_bad, l), std::domain_error,
                   " magnitude");
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

TEST(MathPrimMat, domain_error_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x1(3);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  std::vector<double> x2(4);
  x2[0] = -2;
  x2[1] = -1;
  x2[2] = -0.5;
  x2[3] = -5;

  double sigma_bad = -1;
  double l_bad = -1;

  EXPECT_THROW_MSG(stan::math::cov_exp_quad(x1, x2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::cov_exp_quad(x1, x2, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_2(4);
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

  std::vector<Eigen::Matrix<double, 1, -1> > x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x_rvec_2(4);
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

TEST(MathPrimMat, nan_domain_error_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x1(3);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  std::vector<double> x2(4);
  x2[0] = -2;
  x2[1] = -1;
  x2[2] = -0.5;
  x2[3] = -5;

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();
  double l_bad = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW_MSG(stan::math::cov_exp_quad(x1, x2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::cov_exp_quad(x1, x2, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_2(4);
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

  std::vector<Eigen::Matrix<double, 1, -1> > x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x_rvec_2(4);
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

  std::vector<double> x1_bad(x1);
  x1_bad[1] = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x2_bad(x2);
  x2_bad[1] = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_1_bad(x_vec_1);
  x_vec_1_bad[1](1) = std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_2_bad(x_vec_2);
  x_vec_2_bad[1](1) = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, 1, -1> > x_rvec_1_bad(x_rvec_1);
  x_rvec_1_bad[1](1) = std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, 1, -1> > x_rvec_2_bad(x_rvec_2);
  x_rvec_2_bad[1](1) = std::numeric_limits<double>::quiet_NaN();

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

TEST(MathPrimMat, dim_mismatch_vec_eigen_vec_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(4, 1);
    x_vec_2[i] << 4, 1, 3, 1;
  }
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(MathPrimMat, dim_mismatch_vec_eigen_rvec_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1> > x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(1, 3);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(1, 4);
    x_vec_2[i] << 4, 1, 3, 1;
  }
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(MathPrimMat, dim_mismatch_vec_eigen_mixed_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1> > x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_1(4);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(4, 1);
    x_vec_1[i] << 4, 1, 3, 1;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x_vec_2(3);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x_rvec_2(4);
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
