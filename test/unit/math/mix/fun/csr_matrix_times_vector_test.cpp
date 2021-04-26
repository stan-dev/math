#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, csr_matrix_times_vector) {
  auto f = [](const auto& w, const auto& b) {
    using stan::math::csr_matrix_times_vector;
    std::vector<int> v{1, 2, 0, 2, 4, 2, 1, 4};
    std::vector<int> u{0, 2, 4, 5, 6, 8};
    return csr_matrix_times_vector(5, 5, w, v, u, b);
  };

  Eigen::VectorXd w(8);
  w << 22, 7, 3, 5, 14, 1, 17, 8;

  Eigen::VectorXd b(8);
  b << 1, 2, 3, 4, 5, 6, 7, 8;

  stan::test::expect_ad(f, w, b);
}
