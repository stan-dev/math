#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mixScalFun, conj) {
  auto f = [](const auto& x) {
    using stan::math::conj;
    return conj(x);
  };
  stan::test::expect_complex_common(f);
}

template <typename T>
void test_vectorized_conj() {
  using complex_t = std::complex<T>;
  using matrix_t = Eigen::Matrix<T, -1, -1>;
  using complex_matrix = Eigen::Matrix<complex_t, -1, -1>;
  complex_matrix A(2, 2);
  A << complex_t(T(0), T(1)), complex_t(T(2), T(3)), complex_t(T(4), T(5)),
      complex_t(T(6), T(7));

  auto f = [](const auto& x) {
    using stan::math::conj;
    return conj(x);
  };
  stan::test::expect_ad(f, A);
}

TEST(mixScalFun, conj_vectorized) {
  test_vectorized_conj<double>();
  test_vectorized_conj<stan::math::var>();
}
