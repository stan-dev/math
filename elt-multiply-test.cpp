#include <cmath>
#include <stan/math.hpp>
#include <iostream>


int main() {
	Eigen::Matrix<stan::math::var, 2, 2> M1;
	Eigen::Matrix<stan::math::var, 2, 2> M2;
	Eigen::Matrix<stan::math::var, 2, 2> M_res1;
	M1 << 1,1,1,1;
	M2 << 2,2,2,2;
	M_res1 = elt_multiply(M1, M2);
	stan::math::var lp2 = 0;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			// lp2 += pow(M_var(i,j), 2
			lp2 += M_res1(i,j);
		}
	}
	lp2.grad();
	std::cout << M1.adj() << std::endl;
	
	
	Eigen::Matrix<double, 2, 2> M3;
	Eigen::Matrix<stan::math::var, 2, 2> M4;
	Eigen::Matrix<stan::math::var, 2, 2> M_res2;
	M3 << 1,1,1,1;
	M4 << 2,2,2,2;
	M_res2 = elt_multiply(M3, M4);
	
	
	return 0;
}
