#include <cmath>
#include <stan/math.hpp>
#include <iostream>


int main() {
	Eigen::Matrix<stan::math::var, 1, 2> M1;
	Eigen::Matrix<stan::math::var, 2, 2> M2;
	Eigen::Matrix<stan::math::var, 2, 2> M_res1;
	M1 << 1,1;
	M2 << 2,2,2,2;
	M_res1 = diag_pre_multiply(M1, M2);
	
	
	Eigen::Matrix<double, 1, 2> M3;
	Eigen::Matrix<stan::math::var, 2, 2> M4;
	Eigen::Matrix<stan::math::var, 2, 2> M_res2;
	M3 << 1,1;
	M4 << 2,2,2,2;
	M_res2 = diag_pre_multiply(M3, M4);
	
	
	return 0;
}
