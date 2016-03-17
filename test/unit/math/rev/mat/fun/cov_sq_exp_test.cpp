#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

TEST(RevMath, cov_sq_exp_vvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * sq_distance/(sq_l * l.val()), grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * -distance/sq_l, grad[2])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * distance/sq_l, grad[3])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * -distance/sq_l, grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * distance/sq_l, grad[2])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}


TEST(RevMath, cov_sq_exp_vdv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double sigma = 0.2;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var l = 5;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * sq_distance/(sq_l * l.val()), grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * -distance/sq_l, grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * distance/sq_l, grad[2])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vdd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double sigma = 0.2;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * -distance/sq_l, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * distance/sq_l, grad[1])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_dvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
  
      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * sq_distance/(sq_l * l.val()), grad[1])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_dvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
  
      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_ddv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double sigma = 0.2;
      
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * sq_distance/(sq_l * l.val()), grad[0]);

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vector_vvv) {
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
      
      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * sq_distance/(sq_l * l.val()), grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[i](0).val() - x[j](0).val()) / sq_l,
                      grad[2])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[i](1).val() - x[j](1).val()) / sq_l,
                      grad[3])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[j](0).val() - x[i](0).val()) / sq_l,
                      grad[4])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[j](1).val() - x[i](1).val()) / sq_l,
                      grad[5])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vector_vvd) {
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
      
      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[i](0).val() - x[j](0).val()) / sq_l,
                      grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[i](1).val() - x[j](1).val()) / sq_l,
                      grad[2])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[j](0).val() - x[i](0).val()) / sq_l,
                      grad[3])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[j](1).val() - x[i](1).val()) / sq_l,
                      grad[4])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vector_vdv) {
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
      
      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * sq_distance/(sq_l * l.val()), grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[i](0).val() - x[j](0).val()) / sq_l,
                      grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[i](1).val() - x[j](1).val()) / sq_l,
                      grad[2])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[j](0).val() - x[i](0).val()) / sq_l,
                      grad[3])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[j](1).val() - x[i](1).val()) / sq_l,
                      grad[4])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vector_vdd) {
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

      
      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[i](0).val() - x[j](0).val()) / sq_l,
                      grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[i](1).val() - x[j](1).val()) / sq_l,
                      grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[j](0).val() - x[i](0).val()) / sq_l,
                      grad[2])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(-cov(i, j).val()
                      * (x[j](1).val() - x[i](1).val()) / sq_l,
                      grad[3])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vector_dvv) {
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
      
      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);

      cov(i, j).grad(params, grad);

      // std::cout << "x[0]: " << x[0] << std::endl
      //           << "x[1]: " << x[1] << std::endl
      //           << "x[2]: " << x[2] << std::endl
      //           << std::endl;
      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_distance = stan::math::squared_distance(stan::math::value_of(x[i]),
                                                        stan::math::value_of(x[j]));
      
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val,
                      grad[0])
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

TEST(RevMath, cov_sq_exp_vector_dvd) {
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
      
      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);

      cov(i, j).grad(params, grad);

      double distance = stan::math::distance(stan::math::value_of(x[i]),
                                             stan::math::value_of(x[j]));
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vector_ddv) {
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
      
      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

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
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * sq_distance/(sq_l * l.val()), grad[0])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

