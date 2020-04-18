#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(AgradRevMatrix, multiply_val_vv_cl) {
  stan::math::matrix_v G(3,3); 
  Eigen::VectorXd g0(3);
  Eigen::MatrixXd CE(3,1);
  Eigen::VectorXd ce0(1);
  Eigen::MatrixXd CI(3,4); 
  Eigen::VectorXd ci0(4);
  
  G << 2.1, 1.5, 1.2,
       1.5, 2.2, 1.3,
       1.2, 1.3, 3.1;
  
  g0 << 6, 1, 1;
  
  CE << 1, 2, -1;
  
  ce0(0)=-4;
  
  CI << 1, 0, 0, -1,
        0, 1, 0, -1,
        0, 0, 1,  0;
  
  
  ci0 << 0, 0, 0, 10;
  

  stan::math::vector_v x = solve_quadprog(G, g0,  CE, ce0,  CI, ci0);
}