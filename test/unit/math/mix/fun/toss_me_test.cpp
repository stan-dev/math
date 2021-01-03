#include <test/unit/math/test_ad.hpp>


TEST(Example, ConstCoeffRef) {
  using stan::math::var;
  Eigen::Matrix<var, -1, -1> A = Eigen::MatrixXd::Random(4, 4);
  Eigen::Matrix<var, -1, -1> B = Eigen::MatrixXd::Random(4, 4);
  Eigen::Matrix<double, -1, -1> C = A.val() * B.adj();
}
