#include <cmath>
#include <stan/math.hpp>
#include <iostream>


int main() {
	// Both var ------------------------------------------------------------------
	Eigen::Matrix<stan::math::var, 1, 2> V1;
	Eigen::Matrix<stan::math::var, 2, 1> V2;
	V1 << 1, 2;
	V2 << 3, 4;
	std::cout << V1.replicate(1, 5) << std::endl;
	std::cout << V1.replicate(5, 1) << std::endl;
	std::cout << V1.transpose().replicate(1, 5) << std::endl;
	std::cout << V1.transpose().replicate(5, 1) << std::endl;
	std::cout << V2.replicate(1, 5) << std::endl;
	std::cout << V2.replicate(5, 1) << std::endl;
	return 0;
}
