#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(AgradRevMatrix, dddddd) {
  stan::math::matrix_d G(3,3); 
  stan::math::vector_d g0(3);
  Eigen::MatrixXd CE(3,1);
  Eigen::VectorXd ce0(1);
  Eigen::MatrixXd CI(3,4); 
  Eigen::VectorXd ci0(4);

  G << 2.1, 0.0, 1.0,
        1.5, 2.2, 0.0,
        1.2, 1.3, 3.1;
  
  g0 << 6, 1, 1;
  
  CE << 1, 2, -1;
  
  ce0(0)=-4;
  
  CI << 1, 0, 0, -1,
        0, 1, 0, -1,
        0, 0, 1,  0;
  
  ci0 << 0, 0, 0, 10;


  stan::math::vector_d x = stan::math::solve_quadprog(G, g0,  CE, ce0,  CI, ci0);
  stan::math::vector_d solution(3);
  solution << -3e-16, 2, 0;
  for(int i=0;i<x.size();i++) {
      EXPECT_NEAR(solution(i), x(i), 1e-8);
  }
}

TEST(AgradRevMatrix, vddddd) {
  stan::math::matrix_v G(3,3); 
  stan::math::vector_d g0(3);
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


  stan::math::vector_v x = stan::math::solve_quadprog(G, g0,  CE, ce0,  CI, ci0);
  stan::math::vector_d solution(3);
  solution << -3e-16, 2, 0;
  for(int i=0;i<x.size();i++) {
      EXPECT_NEAR(solution(i), x(i).val(), 1e-8);
  }
}

TEST(AgradRevMatrix, vvdddd) {
  stan::math::matrix_v G(3,3); 
  stan::math::vector_v g0(3);
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


  stan::math::vector_v x = stan::math::solve_quadprog(G, g0,  CE, ce0,  CI, ci0);
  stan::math::vector_d solution(3);
  solution << -3e-16, 2, 0;
  for(int i=0;i<x.size();i++) {
      EXPECT_NEAR(solution(i), x(i).val(), 1e-8);
  }
}

TEST(AgradRevMatrix, vvvddd) {
  stan::math::matrix_v G(3,3); 
  stan::math::vector_v g0(3);
  stan::math::matrix_v CE(3,1);
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


  stan::math::vector_v x = stan::math::solve_quadprog(G, g0,  CE, ce0,  CI, ci0);
  stan::math::vector_d solution(3);
  solution << -3e-16, 2, 0;
  for(int i=0;i<x.size();i++) {
      EXPECT_NEAR(solution(i), x(i).val(), 1e-8);
  }
}

TEST(AgradRevMatrix, vvvvdd) {
  stan::math::matrix_v G(3,3); 
  stan::math::vector_v g0(3);
  stan::math::matrix_v CE(3,1);
  stan::math::vector_v ce0(1);
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


  stan::math::vector_v x = stan::math::solve_quadprog(G, g0,  CE, ce0,  CI, ci0);
  stan::math::vector_d solution(3);
  solution << -3e-16, 2, 0;
  for(int i=0;i<x.size();i++) {
      EXPECT_NEAR(solution(i), x(i).val(), 1e-8);
  }
}

TEST(AgradRevMatrix, vvvvvd) {
  stan::math::matrix_v G(3,3); 
  stan::math::vector_v g0(3);
  stan::math::matrix_v CE(3,1);
  stan::math::vector_v ce0(1);
  stan::math::matrix_v CI(3,4); 
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


  stan::math::vector_v x = stan::math::solve_quadprog(G, g0,  CE, ce0,  CI, ci0);
  stan::math::vector_d solution(3);
  solution << -3e-16, 2, 0;
  for(int i=0;i<x.size();i++) {
      EXPECT_NEAR(solution(i), x(i).val(), 1e-8);
  }
}

TEST(AgradRevMatrix, vvvvvv) {
  stan::math::matrix_v G(3,3); 
  stan::math::vector_v g0(3);
  stan::math::matrix_v CE(3,1);
  stan::math::vector_v ce0(1);
  stan::math::matrix_v CI(3,4); 
  stan::math::vector_v ci0(4);

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


  stan::math::vector_v x = stan::math::solve_quadprog(G, g0,  CE, ce0,  CI, ci0);
  stan::math::vector_d solution(3);
  solution << -3e-16, 2, 0;
  for(int i=0;i<x.size();i++) {
      EXPECT_NEAR(solution(i), x(i).val(), 1e-8);
  }
}