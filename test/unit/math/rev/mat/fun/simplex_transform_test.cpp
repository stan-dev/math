#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/jacobian.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <vector>
#include <random>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(probTransform, simplex_jacobian) {
  using stan::math::var;
  using std::vector;
  var a = 2.0;
  var b = 3.0;
  var c = -1.0;

  Matrix<var, Dynamic, 1> y(3);
  y << a, b, c;

  var lp(0);
  Matrix<var, Dynamic, 1> x = stan::math::simplex_constrain(y, lp);

  vector<var> indeps;
  indeps.push_back(a);
  indeps.push_back(b);
  indeps.push_back(c);

  vector<var> deps;
  deps.push_back(x(0));
  deps.push_back(x(1));
  deps.push_back(x(2));

  vector<vector<double> > jacobian;
  stan::math::jacobian(deps, indeps, jacobian);

  Matrix<double, Dynamic, Dynamic> J(3, 3);
  for (int m = 0; m < 3; ++m) {
    for (int n = 0; n < 3; ++n) {
      J(m, n) = jacobian[m][n];
    }
  }

  double det_J = J.determinant();
  double log_det_J = log(det_J);

  EXPECT_FLOAT_EQ(log_det_J, lp.val());
}

TEST(prob_transform, simplex_constrain_length_zero_no_segfault) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> xv(0);
  var out = sum(simplex_constrain(xv));

  out.grad();
}

TEST(prob_transform, simplex_constrain_length_one_no_segfault) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> xv(1);
  xv << 1.0;
  var out = sum(simplex_constrain(xv));

  out.grad();

  EXPECT_FLOAT_EQ(xv(0).adj(), 0.0);
}

TEST(AgradRevMatrix, simplex_constrain_analytical_gradients) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> yv(4);
  Eigen::Matrix<double, Eigen::Dynamic, 1> zk(yv.size());
  Eigen::Matrix<double, Eigen::Dynamic, 1> diagonal_of_jacobian(yv.size());
  yv << 2.0, -1.0, 0.0, -1.1;
  Eigen::Matrix<var, Eigen::Dynamic, 1> xv = simplex_constrain(yv);
  var out = 0.0;
  for (int i = 0; i < yv.size(); ++i)
    out += xv[i];
  double tk = 1.0;
  int Km1 = yv.size();
  for (int k = 0; k < diagonal_of_jacobian.size(); ++k) {
    zk(k) = stan::math::inv_logit(yv[k].val() - log(Km1 - k));
    diagonal_of_jacobian(k) = tk * zk(k) * (1.0 - zk(k));
    tk -= xv(k).val();
  }

  out.grad();

  EXPECT_FLOAT_EQ(yv(0).adj(), diagonal_of_jacobian(0) * (1 - zk(3))
                                   * (1 - zk(2)) * (1 - zk(1)));
  EXPECT_FLOAT_EQ(yv(1).adj(),
                  diagonal_of_jacobian(1) * (1 - zk(3)) * (1 - zk(2)));
  EXPECT_FLOAT_EQ(yv(2).adj(), diagonal_of_jacobian(2) * (1 - zk(3)));
  EXPECT_FLOAT_EQ(yv(3).adj(), diagonal_of_jacobian(3));

  stan::math::set_zero_all_adjoints();

  out = sum(xv);

  out.grad();

  EXPECT_FLOAT_EQ(yv(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(yv(1).adj(), 0.0);
  EXPECT_FLOAT_EQ(yv(2).adj(), 0.0);
  EXPECT_FLOAT_EQ(yv(3).adj(), 0.0);
}

TEST(prob_transform, simplex_constrain_analytical_grads_rng) {
  using stan::math::var;

  std::mt19937 rng;
  std::uniform_real_distribution<> uniform(-5.0, 5.0);

  for (int i = 0; i < 5; ++i) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> yv(4);
    Eigen::Matrix<double, Eigen::Dynamic, 1> zk(yv.size());
    Eigen::Matrix<double, Eigen::Dynamic, 1> diagonal_of_jacobian(yv.size());
    for (int j = 0; j < yv.size(); j++) {
      yv(j) = uniform(rng);
    }
    Eigen::Matrix<var, Eigen::Dynamic, 1> xv = simplex_constrain(yv);
    var out = 0.0;
    for (int i = 0; i < yv.size(); ++i)
      out += xv[i];
    double tk = 1.0;
    int Km1 = yv.size();
    for (int k = 0; k < diagonal_of_jacobian.size(); ++k) {
      zk(k) = stan::math::inv_logit(yv[k].val() - log(Km1 - k));
      diagonal_of_jacobian(k) = tk * zk(k) * (1.0 - zk(k));
      tk -= xv(k).val();
    }

    out.grad();

    EXPECT_FLOAT_EQ(yv(0).adj(), diagonal_of_jacobian(0) * (1 - zk(3))
                                     * (1 - zk(2)) * (1 - zk(1)));
    EXPECT_FLOAT_EQ(yv(1).adj(),
                    diagonal_of_jacobian(1) * (1 - zk(3)) * (1 - zk(2)));
    EXPECT_FLOAT_EQ(yv(2).adj(), diagonal_of_jacobian(2) * (1 - zk(3)));
    EXPECT_FLOAT_EQ(yv(3).adj(), diagonal_of_jacobian(3));

    stan::math::set_zero_all_adjoints();

    out = sum(xv);

    out.grad();

    EXPECT_FLOAT_EQ(yv(0).adj(), 0.0);
    EXPECT_FLOAT_EQ(yv(1).adj(), 0.0);
    EXPECT_FLOAT_EQ(yv(2).adj(), 0.0);
    EXPECT_FLOAT_EQ(yv(3).adj(), 0.0);
  }
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::vector_v y(3);
  y << 2, 3, -1;
  stan::math::var lp = 0;
  test::check_varis_on_stack(stan::math::simplex_constrain(y, lp));
  test::check_varis_on_stack(stan::math::simplex_constrain(y));
}
