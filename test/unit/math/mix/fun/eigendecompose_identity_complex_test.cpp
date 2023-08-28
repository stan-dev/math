#include <test/unit/math/test_ad.hpp>

template <typename T>
void expect_identity_matrix_complex(const T& x) {
  EXPECT_EQ(x.rows(), x.cols());
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      EXPECT_NEAR(i == j ? 1 : 0, stan::math::value_of_rec(x(i, j).real()),
                  1e-6);
    }
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
            * eigenvalues.asDiagonal().inverse());

  expect_identity_matrix_complex(I);

  std::tie(eigenvectors, eigenvalues) = stan::math::eigendecompose(c22);
  auto I2 = (eigenvectors.inverse() * c22 * eigenvectors
             * eigenvalues.asDiagonal().inverse());

  expect_identity_matrix_complex(I2);
}

TEST(mathMixFun, eigenvectorsIdComplex) {
  using d_t = double;
  using v_t = stan::math::var;
  using fd_t = stan::math::fvar<double>;
  using ffd_t = stan::math::fvar<fd_t>;
  using fv_t = stan::math::fvar<stan::math::var>;
  using ffv_t = stan::math::fvar<fv_t>;

  expectComplexEigenvectorsId<d_t>();
  expectComplexEigenvectorsId<v_t>();
  expectComplexEigenvectorsId<fd_t>();
  expectComplexEigenvectorsId<ffd_t>();
  expectComplexEigenvectorsId<fv_t>();
  expectComplexEigenvectorsId<ffv_t>();
}
