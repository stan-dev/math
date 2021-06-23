#include <test/unit/math/test_ad.hpp>

TEST(mathMixFun, eigenvalues2) {
  auto f = [](const auto& x) {
    using stan::math::eigenvalues;
    return eigenvalues(x);
  };

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
}

// see eigenvectors_test.cpp for test of eigenvectors() and eigenvalues()
// using reconstruction identities
