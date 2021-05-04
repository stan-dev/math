#include <test/unit/math/test_ad.hpp>
#include <stdexcept>

TEST(mathMixFun, pseudoEigenvectors) {
  auto f = [](const auto& x) {
    using stan::math::pseudo_eigenvectors;
    return pseudo_eigenvectors(x);
  };
  for (const auto& x : stan::test::square_test_matrices(0, 2)) {
    stan::test::expect_ad(f, x);
  }

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
}

template <typename T>
void expect_zero_matrix(const T& m) {
  for (int j = 0; j < m.cols(); ++j) {
    for (int i = 0; i < m.rows(); ++i) {
      EXPECT_NEAR(0.0, stan::math::value_of_rec(m(i, j)), 1e-6);
    }
  }
}

template <typename T>
void test_pseudo_eigendecomposition() {
  using stan::math::pseudo_eigenvalues;
  using stan::math::pseudo_eigenvectors;
  for (const auto& x : stan::test::square_test_matrices(1, 3)) {
    Eigen::Matrix<T, -1, -1> A(x);
    auto D = pseudo_eigenvalues(A);
    auto V = pseudo_eigenvectors(A);
    expect_zero_matrix((A - V * D * V.inverse()).eval());
  }
}

TEST(mathMixFun, pseudoEigenVectors) {
  using d_t = double;
  using v_t = stan::math::var;
  using fd_t = stan::math::fvar<d_t>;
  using ffd_t = stan::math::fvar<fd_t>;
  using fv_t = stan::math::fvar<v_t>;
  using ffv_t = stan::math::fvar<fv_t>;
  test_pseudo_eigendecomposition<d_t>();
  test_pseudo_eigendecomposition<v_t>();
  test_pseudo_eigendecomposition<fd_t>();
  test_pseudo_eigendecomposition<ffd_t>();
  test_pseudo_eigendecomposition<fv_t>();
  test_pseudo_eigendecomposition<ffv_t>();
}
