#include <cmath>
#include <stan/math.hpp>
#include <iostream>


int main() {
	// Both var ------------------------------------------------------------------
	Eigen::Matrix<stan::math::var, 1, 2> M1;
	Eigen::Matrix<stan::math::var, 2, 2> M2;
	Eigen::Matrix<stan::math::var, 2, 2> M_res1;
	M1 << 1,1;
	M2 << 2,2,2,2;
	M_res1 = diag_pre_multiply(M1, M2);
	std::cout << M1.adj() << std::endl;
	
	stan::math::var lp;
	// lp = 0; // Ce ta vrstica ni zakomentirana, potem crashne.

	Eigen::Matrix<stan::math::var, 1, 2> M3;
	Eigen::Matrix<stan::math::var, 2, 2> M4;
	Eigen::Matrix<stan::math::var, 2, 2> M_res2;
	M3 << 1,1;
	M4 << 2,2,2,2;
	M_res2 = diag_pre_multiply(M3, M4);
	
	stan::math::var lp2 = 0;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			// lp2 += pow(M_var(i,j), 2);
			lp2 += M_res2(i,j);
		}
	}
	lp2.grad();
	std::cout << M3.adj() << std::endl;

	return 0;
}
