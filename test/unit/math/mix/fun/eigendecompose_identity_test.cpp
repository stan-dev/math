#include <test/unit/math/test_ad.hpp>

template <typename T>
void expect_identity_matrix(const T& x) {
  EXPECT_EQ(x.rows(), x.cols());
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      EXPECT_NEAR(i == j ? 1 : 0, stan::math::value_of_rec(x(i, j)), 1e-6);
    }
  }
}

template <typename T>
void expectEigenvectorsId() {
  for (const auto& m_d : stan::test::square_test_matrices(1, 2)) {
    Eigen::Matrix<T, -1, -1> m(m_d);
    auto vecs = eigenvectors(m).eval();
    auto vals = eigenvalues(m).eval();
    auto I = (vecs.inverse() * m * vecs * vals.asDiagonal().inverse()).real();
    expect_identity_matrix(I);
  }
}

TEST(mathMixFun, eigenvectorsId) {
  using d_t = double;
  using v_t = stan::math::var;
  using fd_t = stan::math::fvar<double>;
  using ffd_t = stan::math::fvar<fd_t>;
  using fv_t = stan::math::fvar<stan::math::var>;
  using ffv_t = stan::math::fvar<fv_t>;

  expectEigenvectorsId<v_t>();
  expectEigenvectorsId<fd_t>();
  expectEigenvectorsId<ffd_t>();
  expectEigenvectorsId<fv_t>();
  expectEigenvectorsId<ffv_t>();
}
