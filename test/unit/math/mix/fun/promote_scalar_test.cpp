#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <complex>
#include <vector>
#include <type_traits>

namespace stan {
namespace math {
namespace test {
template <typename T, typename S, require_all_stan_scalar_t<T, S>* = nullptr>
void expect_same_value_of_rec(T&& a, S&& b) {
  using stan::math::value_of_rec;
  EXPECT_EQ(value_of_rec(a), value_of_rec(b));
}

template <typename T, typename S, require_all_eigen_t<T, S>* = nullptr>
void expect_same_value_of_rec(T&& a, S&& b) {
  using stan::math::value_of_rec;
  EXPECT_MATRIX_EQ(value_of_rec(a), value_of_rec(b));
}

template <typename T, typename S, require_all_tuple_t<T, S>* = nullptr>
void expect_same_value_of_rec(T&& a, S&& b);

template <typename T, typename S, require_all_std_vector_t<T, S>* = nullptr>
void expect_same_value_of_rec(T&& a, S&& b) {
  using stan::math::value_of_rec;
  for (size_t i = 0; i < a.size(); ++i) {
    expect_same_value_of_rec(a[i], b[i]);
  }
}

template <typename T, typename S, require_all_tuple_t<T, S>*>
void expect_same_value_of_rec(T&& a, S&& b) {
  stan::math::for_each(
      [](auto&& test, auto&& result) {
        stan::math::test::expect_same_value_of_rec(test, result);
      },
      a, b);
}

}  // namespace test
}  // namespace math
}  // namespace stan

template <typename PromotionType, typename UnPromotedType>
void test_promote_scalar() {
  using eig_mat = Eigen::Matrix<UnPromotedType, -1, -1>;
  auto tester_gen = []() {
    return std::make_tuple(UnPromotedType(1.0), eig_mat::Random(2, 2).eval(),
                           std::vector<UnPromotedType>{1, 2, 3},
                           std::vector<eig_mat>{eig_mat::Random(2, 2).eval(),
                                                eig_mat::Random(2, 2).eval()});
  };
  using tester_t
      = std::tuple<UnPromotedType, eig_mat, std::vector<UnPromotedType>,
                   std::vector<eig_mat>>;
  auto tester = std::tuple_cat(
      tester_gen(),
      std::make_tuple(std::vector<tester_t>{tester_gen(), tester_gen()}));
  using inner_promo_type
      = std::tuple<PromotionType, PromotionType, PromotionType, PromotionType>;
  auto result = stan::math::promote_scalar<
      std::tuple<PromotionType, PromotionType, PromotionType, PromotionType,
                 inner_promo_type>>(tester);
  stan::math::test::expect_same_value_of_rec(tester, result);
}

TEST(mixFun, promote_scalar_tuple) {
  using stan::math::fvar;
  using stan::math::var;
  test_promote_scalar<double, int>();
  test_promote_scalar<double, double>();
  test_promote_scalar<var, int>();
  test_promote_scalar<var, double>();
  test_promote_scalar<var, var>();
  test_promote_scalar<fvar<double>, int>();
  test_promote_scalar<fvar<double>, double>();
  test_promote_scalar<fvar<double>, fvar<double>>();
  test_promote_scalar<fvar<var>, double>();
  test_promote_scalar<fvar<var>, var>();
  test_promote_scalar<fvar<var>, fvar<var>>();
  test_promote_scalar<std::complex<double>, double>();
  test_promote_scalar<std::complex<var>, double>();
  test_promote_scalar<std::complex<fvar<var>>, double>();
}

template <typename UnPromotedType>
void test_promote_scalar_basic() {
  std::tuple<UnPromotedType, UnPromotedType> x{UnPromotedType(3.5),
                                               UnPromotedType(4.5)};
  std::tuple<std::complex<UnPromotedType>, UnPromotedType> z
      = stan::math::promote_scalar<
          std::tuple<std::complex<UnPromotedType>, UnPromotedType>>(x);
  stan::math::test::expect_same_value_of_rec(x, z);
}

TEST(mixFun, promote_scalar_tuple_basic) {
  using stan::math::fvar;
  using stan::math::var;
  test_promote_scalar_basic<double>();
  test_promote_scalar_basic<var>();
  test_promote_scalar_basic<fvar<double>>();
  test_promote_scalar_basic<fvar<var>>();
}
