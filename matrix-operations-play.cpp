#include <cmath>
#include <stan/math.hpp>
#include <iostream>


int main() {
	// Both var ------------------------------------------------------------------
	Eigen::Matrix<stan::math::var, 1, 2> M1;
	Eigen::Matrix<stan::math::var, 2, 2> M2;
	Eigen::Matrix<stan::math::var, 2, 2> M_res1;
  Eigen::Matrix<stan::math::var, 2, 1> M4;
	M1 << 1,1;
	M2 << 2,2,2,2;
	M_res1 = diag_pre_multiply(M1, M2);
	std::cout << M1 << std::endl;
  std::cout << M1.val() << std::endl;
  std::cout << M1.val().array() << std::endl;
  M4 = M2.rowwise().sum();
  std::cout << M4 << std::endl;
  // std::cout << M2.rowwise().sum << std::endl;
  std::cout << M2.val().cwiseProduct(M2.val()) << std::endl;
  std::cout << M2.val().cwiseProduct(M2.adj()) << std::endl;
  // std::cout << M1.coeff() << std::endl;
  Eigen::Matrix<stan::math::var, 2, 1> M5;
  M5 =  M2.val().cwiseProduct(M2.val()).rowwise().sum();
  std::cout << M5 << std::endl;

  Eigen::Matrix<stan::math::var, 2, 1> M6;
  M6 =  M2.val().cwiseProduct(M2.adj()).rowwise().sum();
  std::cout << M6 << std::endl;

  std::cout << "Replicate: " << std::endl;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> M_rep;
  M_rep = M1.replicate(8,1);
  std::cout << M_rep << std::endl;
	return 0;
}
