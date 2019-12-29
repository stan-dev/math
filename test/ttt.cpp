#include <stan/math/mix/mat.hpp>
#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/scal/prob/von_mises_lpdf.hpp>
#include <iostream>
/*
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>
#include <string>
*/

int main() {
  double y = -0.8;
  double mu = 0.4;
  double kappa = 0;

  stan::math::von_mises_lpdf(y, mu, kappa);

  /*
    const double t = 1.0;
    Eigen::MatrixXd A(0, 0);
    Eigen::MatrixXd B(0, 0);

    Eigen::Matrix<var, -1, -1> C(0, 2);
    auto D = stan::math::scale_matrix_exp_multiply(t, A, C);
    std::cout << "rows: " << D.rows() << ", cols: " << D.cols() << std::endl;
    double nan = std::numeric_limits<double>::quiet_NaN();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y(1, 1);
    y << nan;

    stan::math::check_pos_definite("ciao", "y", y);
  */
}
