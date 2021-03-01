#include <cmath>
#include <stan/math.hpp>
#include <iostream>


int main() {
	Eigen::Matrix<double, 1, 2> V1;
	Eigen::Matrix<double, 2, 2> M2;
	V1 << 2, 3;
	M2 << 3, 3, 1, 1;
	std::cout << V1 << std::endl << std::endl;
	std::cout << M2 << std::endl << std::endl;
	Eigen::Matrix<double, 2, 2> M3;
	M3 = stan::math::diag_pre_multiply(V1, M2);
	std::cout << M3 << std::endl;
	Eigen::Matrix<double, 2, 2> M4;
	M4 = V1.asDiagonal();
	std::cout << M4 << std::endl;
	Eigen::Matrix<double, 2, 2> M5;
	M5 = V1.asDiagonal() * M2;
	std::cout << M5 << std::endl;
	
	
	// First one is var ----------------------------------------------------------
	Eigen::Matrix<stan::math::var, 1, 2> V2;
	Eigen::Matrix<stan::math::var, 2, 2> M_var;
	V2 << 2, 3;
	// stan::math::diag_pre_multiply(V2, M2); // For now, since it is void.
	/* M_var = stan::math::diag_pre_multiply(V2, M2);
	
	stan::math::var lp2 = 0;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			// lp2 += pow(M_var(i,j), 2
			lp2 += M_var(i,j);
		}
	}
	lp2.grad();
	std::cout << V2.adj() << std::endl; // Expect 72, 12 */
	
	// Second one is var ---------------------------------------------------------
	/* Eigen::Matrix<double, 1, 2> V3;
	Eigen::Matrix<stan::math::var, 2, 2> M6;
	Eigen::Matrix<stan::math::var, 2, 2> M_var2;
	V3 << 2, 1;
	M6 << 2, 3, 1, 4;
	M_var2 = stan::math::diag_pre_multiply(V3, M6);
	
	stan::math::var lp3 = 0;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			// lp2 += pow(M_var(i,j), 2
			lp3 += M_var2(i,j);
		}
	}
	lp3.grad();
	std::cout << M6.adj() << std::endl; */
	
	// Both are var --------------------------------------------------------------
	Eigen::Matrix<stan::math::var, 1, 2> V_both;
	Eigen::Matrix<stan::math::var, 2, 2> M_both;
	Eigen::Matrix<stan::math::var, 2, 2> M_both_res;
	V_both << 2, 4;
	M_both << 1, 2, 3, 4;
	M_both_res = stan::math::diag_pre_multiply(V_both, M_both);
	stan::math::var lp4 = 0;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			// lp2 += pow(M_var(i,j), 2
			lp4 += M_both_res(i,j);
		}
	}
	lp4.grad();
	std::cout << M_both.adj() << std::endl << std::endl;
	std::cout << V_both.adj() << std::endl;
	
	return 0;
}
