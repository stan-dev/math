#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathPrim, basic_print) {
  int i = 1;
  double d = 1.0;
  Eigen::VectorXd v = Eigen::VectorXd::Ones(1);
  Eigen::RowVectorXd rv = Eigen::RowVectorXd::Ones(1);
  Eigen::MatrixXd m = Eigen::MatrixXd::Ones(1, 1);
  std::vector<double> a(1, 1);

  std::tuple<Eigen::VectorXd, int, std::vector<double>> tup(v, i, a);

  {
    std::stringstream s;
    stan::math::stan_print(&s, i);
    EXPECT_TRUE(s.str().find("1") != std::string::npos);
  }

  {
    std::stringstream s;
    stan::math::stan_print(&s, d);
    EXPECT_TRUE(s.str().find("1") != std::string::npos);
  }

  {
    std::stringstream s;
    stan::math::stan_print(&s, v);
    EXPECT_TRUE(s.str().find("[1]") != std::string::npos);
  }

  {
    std::stringstream s;
    stan::math::stan_print(&s, rv);
    EXPECT_TRUE(s.str().find("[1]") != std::string::npos);
  }

  {
    std::stringstream s;
    stan::math::stan_print(&s, m);
    EXPECT_TRUE(s.str().find("[[1]]") != std::string::npos);
  }

  {
    std::stringstream s;
    stan::math::stan_print(&s, a);
    EXPECT_TRUE(s.str().find("[1]") != std::string::npos);
  }

  {
    std::stringstream s;
    stan::math::stan_print(&s, tup);
    EXPECT_TRUE(s.str().find("([1],1,[1])") != std::string::npos);
  }
}

TEST(MathPrim, basic_expressions) {
  Eigen::VectorXd v = Eigen::VectorXd::Ones(2);
  Eigen::RowVectorXd rv = Eigen::RowVectorXd::Ones(2);
  Eigen::MatrixXd m = Eigen::MatrixXd::Ones(2, 2);

  {
    std::stringstream s;
    stan::math::stan_print(&s, v * v.transpose());
    EXPECT_TRUE(s.str().find("[[1,1],[1,1]]") != std::string::npos);
  }

  {
    std::stringstream s;
    stan::math::stan_print(&s, rv * m);
    EXPECT_TRUE(s.str().find("[2,2]") != std::string::npos);
  }

  {
    std::stringstream s;
    stan::math::stan_print(&s, m * m);
    EXPECT_TRUE(s.str().find("[[2,2],[2,2]]") != std::string::npos);
  }
}
