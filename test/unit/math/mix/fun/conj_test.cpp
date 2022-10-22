#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mathMixMatFun, conj) {
  auto f = [](const auto& x) { return stan::math::conj(x); };
  stan::test::expect_complex_common(f);
  stan::test::expect_unary_vectorized<stan::test::ScalarSupport::ComplexOnly>(
      f);
}

template <typename T>
void test_vectorized_conj() {
  using stan::math::value_of_rec;
  using complex_t = std::complex<T>;
  using complex_matrix = Eigen::Matrix<complex_t, -1, -1>;
  complex_matrix A(2, 2);
  A << complex_t(T(0), T(1)), complex_t(T(2), T(3)), complex_t(T(4), T(5)),
      complex_t(T(6), T(7));
  auto A_conj = stan::math::conj(A);
  EXPECT_MATRIX_COMPLEX_FLOAT_EQ(value_of_rec(A_conj),
                                 value_of_rec(A.conjugate()))
  std::vector<complex_t> v{complex_t(T(0), T(1)), complex_t(T(2), T(3)),
                           complex_t(T(4), T(5)), complex_t(T(6), T(7))};

  std::vector<complex_t> v_conj = stan::math::conj(v);
  for (int i = 0; i < v.size(); ++i) {
    std::complex<double> vi = value_of_rec(v_conj[i]);
    std::complex<double> ci = value_of_rec(stan::math::conj(v[i]));
    EXPECT_FLOAT_EQ(vi.real(), ci.real());
    EXPECT_FLOAT_EQ(vi.imag(), ci.imag());
  }
}

TEST(mathMixMatFun, conj_vectorized) {
  using d_t = double;
  using v_t = stan::math::var;
  using fd_t = stan::math::fvar<d_t>;
  using ffd_t = stan::math::fvar<fd_t>;
  using fv_t = stan::math::fvar<v_t>;
  using ffv_t = stan::math::fvar<fv_t>;

  test_vectorized_conj<d_t>();
  test_vectorized_conj<v_t>();
  test_vectorized_conj<fd_t>();
  test_vectorized_conj<ffd_t>();
  test_vectorized_conj<fv_t>();
  test_vectorized_conj<ffv_t>();
}
