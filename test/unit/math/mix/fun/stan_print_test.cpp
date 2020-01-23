#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

TEST(AgradMixMatrixStanPrint, fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::stan_print;
  using stan::math::var;
  using std::vector;

  std::stringstream output;
  fvar<var> a(1, 2);
  stan_print(&output, a);
  EXPECT_EQ("1", output.str());

  output.str(std::string());
  output.clear();
  std::vector<fvar<var> > b;
  b.push_back(a);
  b.push_back(a);
  b.push_back(a);
  stan_print(&output, b);
  EXPECT_EQ("[1,1,1]", output.str());

  output.str(std::string());
  output.clear();
  Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> c(3);
  c << a, a, a;
  stan_print(&output, c);
  EXPECT_EQ("[1,1,1]", output.str());

  output.str(std::string());
  output.clear();
  Eigen::Matrix<fvar<var>, 1, Eigen::Dynamic> d(3);
  d << a, a, a;
  stan_print(&output, d);
  EXPECT_EQ("[1,1,1]", output.str());

  output.str(std::string());
  output.clear();
  Eigen::Matrix<fvar<var>, Eigen::Dynamic, Eigen::Dynamic> e(2, 2);
  e << a, a, a, a;
  stan_print(&output, e);
  EXPECT_EQ("[[1,1],[1,1]]", output.str());
}

TEST(AgradMixMatrixStanPrint, fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::stan_print;
  using stan::math::var;
  using std::vector;

  std::stringstream output;
  fvar<fvar<var> > a(1, 2);
  stan_print(&output, a);
  EXPECT_EQ("1", output.str());

  output.str(std::string());
  output.clear();
  std::vector<fvar<fvar<var> > > b;
  b.push_back(a);
  b.push_back(a);
  b.push_back(a);
  stan_print(&output, b);
  EXPECT_EQ("[1,1,1]", output.str());

  output.str(std::string());
  output.clear();
  Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1> c(3);
  c << a, a, a;
  stan_print(&output, c);
  EXPECT_EQ("[1,1,1]", output.str());

  output.str(std::string());
  output.clear();
  Eigen::Matrix<fvar<fvar<var> >, 1, Eigen::Dynamic> d(3);
  d << a, a, a;
  stan_print(&output, d);
  EXPECT_EQ("[1,1,1]", output.str());

  output.str(std::string());
  output.clear();
  Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, Eigen::Dynamic> e(2, 2);
  e << a, a, a, a;
  stan_print(&output, e);
  EXPECT_EQ("[[1,1],[1,1]]", output.str());
}
