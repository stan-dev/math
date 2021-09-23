#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, csr_matrix_times_vector) {
  auto f = [](const auto& w, const auto& b) {
    using stan::math::csr_matrix_times_vector;
    std::vector<int> v{1, 2, 3, 1, 2};
    std::vector<int> u{1, 2, 3, 4, 5, 6};
    return csr_matrix_times_vector(5, 5, w, v, u, b);
  };

  Eigen::VectorXd w(5);
  w << -0.67082, 0.5, -0.223607, -0.223607, -0.5;
  Eigen::VectorXd b(5);
  b << 1, 2, 3, 4, 5;

  stan::test::expect_ad(f, w, b);
  stan::test::expect_ad_matvar(f, w, b);
}
