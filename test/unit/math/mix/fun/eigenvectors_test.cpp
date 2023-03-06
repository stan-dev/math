#include <test/unit/math/test_ad.hpp>
#include <stdexcept>

TEST(mathMixFun, eigenvectors) {
  auto f = [](const auto& x) {
    using stan::math::eigenvectors;
    return eigenvectors(x);
  };
  for (const auto& x : stan::test::square_test_matrices(0, 2)) {
    stan::test::expect_ad(f, x);
  }

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
}

TEST(mathMixFun, eigenvectorsComplex) {
  auto f = [](const auto& x) {
    using stan::math::eigenvectors;
    return eigenvectors(stan::math::to_complex(x, 0));
  };
  for (const auto& x : stan::test::square_test_matrices(0, 2)) {
    stan::test::expect_ad(f, x);
  }

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
}

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

template <typename T>
void expectComplexEigenvectorsId() {
  Eigen::Matrix<std::complex<T>, -1, -1> c22(2, 2);
  c22 << stan::math::to_complex(T(0), T(-1)),
      stan::math::to_complex(T(0), T(0)), stan::math::to_complex(T(2), T(0)),
      stan::math::to_complex(T(4), T(0));
  auto eigenvalues = stan::math::eigenvalues(c22);
  auto eigenvectors = stan::math::eigenvectors(c22);

  auto I = (eigenvectors.inverse() * c22 * eigenvectors
            * eigenvalues.asDiagonal().inverse())
               .real();

  expect_identity_matrix(I);
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

  expectComplexEigenvectorsId<v_t>();
  expectComplexEigenvectorsId<fd_t>();
  expectComplexEigenvectorsId<ffd_t>();
  expectComplexEigenvectorsId<fv_t>();
  expectComplexEigenvectorsId<ffv_t>();
}
