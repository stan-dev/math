#include <cmath>
#include <stan/math.hpp>
#include <iostream>
#include <chrono>


int main() {
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;

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
	std::cout << "------------" << std::endl;
	std::cout << V2.replicate(5, 1) << std::endl;
	std::cout << V1.transpose().
								replicate(5, 1) << std::endl;

	Eigen::MatrixXd V_large = Eigen::MatrixXd::Random(1, 10000);
	Eigen::MatrixXd M_large = Eigen::MatrixXd::Random(10000, 10000);
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M_large_res1;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M_large_res2;
	// M_large_res1 = V_large.transpose().replicate(1, 10000).cwiseProduct(M_large);

	auto t1 = high_resolution_clock::now();
  M_large_res1 = V_large.transpose().replicate(1, 10000).cwiseProduct(M_large);
  auto t2 = high_resolution_clock::now();

	auto ms_int = duration_cast<milliseconds>(t2 - t1);
  std::cout << ms_int.count() << "ms\n";


	auto t3 = high_resolution_clock::now();
	for (int i = 0; i < M_large.rows(); ++i) {
		for (int j = 0; j < M_large.cols(); ++j) {
			M_large_res1(i,j) = V_large(i) * M_large(i,j);
		}
	}
	auto t4 = high_resolution_clock::now();
	auto ms_int2 = duration_cast<milliseconds>(t4 - t3);
	std::cout << ms_int2.count() << "ms\n";

	return 0;
}
