#include <iostream>
#include <lib/eigen_3.3.9/Eigen/Dense>
// #include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/symmetrize_from_lower_tri.hpp>
#include <stan/math/prim/fun/symmetrize_from_upper_tri.hpp>

 
using Eigen::MatrixXd;
using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::symmetrize_from_lower_tri;
using stan::math::symmetrize_from_upper_tri;
using std::cout;
 
int main()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;
	Matrix<double, Dynamic, Dynamic> m2(2, 2);
	m2(0,0) = 1;
	m2(1,0) = 2;
	m2(1,1) = 1;
	cout << m2 << std::endl;
	std::cout << symmetrize_from_lower_tri(m2) << std::endl;
	
	Matrix<double, Dynamic, Dynamic> m3(2, 2);
	m3(0,0) = 1;
	m3(0,1) = 2;
	m3(1,1) = 1;
	m3(1,0) = 0;
	cout << m3 << std::endl;
	std::cout << symmetrize_from_upper_tri(m3) << std::endl;
	// Matrix<double, Dynamic, Dynamic> m3(2, 2);
	// m3(2,2) = stan::math::symmetrize_from_lower_tri(m2);
	// std::cout << m3 << std::endl;
}
