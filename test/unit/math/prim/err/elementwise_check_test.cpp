#include <test/unit/util.hpp>
#include <stan/math/prim.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

namespace elementwise_check_test {
auto p = [](const auto& x) { return !stan::math::is_nan(x); };

template <typename T>
void do_check(const T& x) {
  stan::math::elementwise_check([](const auto& x) { return p(x); },
                                "elementwise_check_tests", "x", x, "not nan");
}

template <typename T>
bool do_is(const T& x) {
  return stan::math::elementwise_is([](const auto& x) { return p(x); }, x);
}

template <typename T>
void expect_good(const T& x) {
  EXPECT_NO_THROW(do_check(x));
  EXPECT_TRUE(do_is(x));
}

template <typename T>
void expect_bad(const T& x) {
  EXPECT_THROW(do_check(x), std::domain_error);
  EXPECT_FALSE(do_is(x));
}
}  // namespace elementwise_check_test

TEST(elementwise_check, checks_scalars) {
  elementwise_check_test::expect_good(0);
  elementwise_check_test::expect_bad(stan::math::NOT_A_NUMBER);
}

TEST(elementwise_check, works_elementwise_on_arrays) {
  const double good = 0;
  const double bad = stan::math::NOT_A_NUMBER;
  using v = std::vector<double>;
  using vv = std::vector<v>;
  using vvv = std::vector<vv>;
  elementwise_check_test::expect_good(v{});
  elementwise_check_test::expect_good(vv{});
  elementwise_check_test::expect_good(vvv{});
  elementwise_check_test::expect_good(v{good});
  elementwise_check_test::expect_good(v{good, good, good});
  elementwise_check_test::expect_good(
      vv{v{good, good}, v{good, good}, v{}, v{good}, v{good, good, good}});
  elementwise_check_test::expect_good(
      vvv{vv{v{good, good, good}, v{good, good}, v{}},
          vv{v{good, good}, v{good, good, good, good}, v{good, good, good}},
          vv{v{}, v{good}, v{good, good, good}}, vv{}, vv{}});
  elementwise_check_test::expect_bad(v{bad});
  elementwise_check_test::expect_bad(v{good, bad, good});
  elementwise_check_test::expect_bad(v{bad, bad, bad});
  elementwise_check_test::expect_bad(
      vv{v{good, good}, v{good, good}, v{}, v{bad}, v{good, good, good}});
  elementwise_check_test::expect_bad(
      vvv{vv{v{good, good, good}, v{good, good, good}, v{}},
          vv{v{good, good}, v{good, good}, v{good, good, good}},
          vv{v{}, v{good}, v{good, good, good, good, bad}}, vv{}, vv{}});
}

TEST(elementwise_check, works_on_eigen_types) {
  const double bad = stan::math::NOT_A_NUMBER;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  elementwise_check_test::expect_good(VectorXd{0});
  elementwise_check_test::expect_good(RowVectorXd{0});
  elementwise_check_test::expect_good(MatrixXd{3, 0});
  elementwise_check_test::expect_good(MatrixXd{0, 3});
  elementwise_check_test::expect_good(MatrixXd{0, 0});
  for (int r = 1; r <= 3; ++r) {
    for (int i = 0; i < r; ++i) {
      const VectorXd good_ev = VectorXd::Zero(r);
      VectorXd bad_ev{good_ev};
      bad_ev[i] = bad;
      elementwise_check_test::expect_good(good_ev);
      elementwise_check_test::expect_bad(bad_ev);
      elementwise_check_test::expect_good(Eigen::RowVectorXd{good_ev});
      elementwise_check_test::expect_bad(Eigen::RowVectorXd{bad_ev});
      for (int c = 1; c <= 3; ++c) {
        for (int j = 0; j < c; ++j) {
          const MatrixXd good_m = MatrixXd::Zero(r, c);
          MatrixXd bad_m{good_m};
          bad_m(i, j) = bad;
          elementwise_check_test::expect_good(good_m);
          elementwise_check_test::expect_bad(bad_m);
        }
      }
    }
  }
}

TEST(elementwise_check, works_on_weird_eigen_types) {
  const double bad = stan::math::NOT_A_NUMBER;
  // Static size and expression templates.
  elementwise_check_test::expect_good(Eigen::Matrix<double, 3, 3>::Zero());
  elementwise_check_test::expect_bad(
      Eigen::Matrix<double, 3, 3>::Constant(bad));
  elementwise_check_test::expect_good(Eigen::VectorXd::Zero(3)
                                      + Eigen::VectorXd::Zero(3));
  elementwise_check_test::expect_bad(
      Eigen::Matrix<double, 3, 3>::Constant(bad)
      + Eigen::Matrix<double, 3, 3>::Constant(bad));
  elementwise_check_test::expect_good(Eigen::Vector3d::Zero()
                                      + Eigen::Vector3d::Zero());
  elementwise_check_test::expect_bad(Eigen::Vector3d::Zero()
                                     + Eigen::Vector3d::Constant(bad));
}

TEST(elementwise_check, works_on_a_ragged_mess_of_dynamic_matrices) {
  const double bad = stan::math::NOT_A_NUMBER;
  using m = Eigen::MatrixXd;
  using v = std::vector<m>;
  using vv = std::vector<v>;
  using vvv = std::vector<vv>;
  const m good00{0, 0};
  const m good31 = Eigen::MatrixXd::Zero(3, 1);
  const m good13 = Eigen::MatrixXd::Zero(1, 3);
  m bad31{good31};
  bad31(1, 0) = bad;
  m bad13{good13};
  bad13(0, 1) = bad;
  elementwise_check_test::expect_good(
      vvv{vv{v{}}, vv{v{}, v{good00, good31}}, vv{v{good31}, v{good31, good31}},
          vv{}, vv{v{good00}, v{good00}, v{good13}}, vv{v{good13}}});
  elementwise_check_test::expect_bad(
      vvv{vv{v{}}, vv{v{}, v{good00, good31}}, vv{v{bad31}, v{good31, good31}},
          vv{}, vv{v{good00}, v{good00}, v{good13}}, vv{v{good13}}});
  elementwise_check_test::expect_bad(
      vvv{vv{v{}}, vv{v{}, v{good00, good31}}, vv{v{good31}, v{good31, good31}},
          vv{}, vv{v{good00}, v{good00}, v{good13}}, vv{v{bad13}}});
}

TEST(elementwise_check, error_messages_look_good) {
  const double bad = stan::math::NOT_A_NUMBER;
  Eigen::VectorXd bad_eigen_v = Eigen::VectorXd::Zero(3);
  bad_eigen_v[1] = bad;
  EXPECT_THROW_MSG(elementwise_check_test::do_check(bad_eigen_v),
                   std::domain_error, "[2]");
  Eigen::MatrixXd bad_m = Eigen::MatrixXd::Zero(3, 4);
  bad_m(1, 2) = bad;
  EXPECT_THROW_MSG(elementwise_check_test::do_check(bad_m), std::domain_error,
                   "[2, 3]");
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajorBit>
      bad_mr = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                             Eigen::RowMajorBit>::Zero(3, 4);
  bad_mr(1, 2) = bad;
  EXPECT_THROW_MSG(elementwise_check_test::do_check(bad_mr), std::domain_error,
                   "[2, 3]");
  std::vector<std::vector<double> > bad_vv{std::vector<double>{},
                                           std::vector<double>{bad}};
  EXPECT_THROW_MSG(elementwise_check_test::do_check(bad_vv), std::domain_error,
                   "[2][1]");
}
