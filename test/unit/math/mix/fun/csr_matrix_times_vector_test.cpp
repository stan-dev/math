#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, csr_matrix_times_vector) {
  auto f = [](const auto& w, const auto& b) {
    using stan::math::csr_matrix_times_vector;
    std::vector<int> v{1,2,3,1,2};
    std::vector<int> u{1,2,3,4,5,6};
    return csr_matrix_times_vector(5, 5, w, v, u, b);
  };

  Eigen::VectorXd w(5);
  w << -0.67082, 0.5, -0.223607, -0.223607, -0.5;

  Eigen::VectorXd b(5);
  b << 1, 2, 3, 4, 5;

  stan::test::expect_ad(f, w, b);
  stan::test::expect_ad_matvar(f, w, b);
}

TEST(MathMixMatFun, csr_matrix_times_vector_err) {

  Eigen::VectorXd w(8);
  w << 22, 7, 3, 5, 14, 1, 17, 8;

  Eigen::VectorXd b(8);
  b << 1, 2, 3, 4, 5, 6, 7, 8;
  std::vector<int> v_too_long{2, 3, 1, 3, 5, 3, 2, 5};
  std::vector<int> u{1, 3, 5, 7, 7, 9};
  EXPECT_THROW(stan::math::csr_matrix_times_vector(8, 8, w, v_too_long, u, b), std::invalid_argument);
  std::vector<int> v_zero{0, 3, 4, 2, 8, 3, 4};
  EXPECT_THROW(stan::math::csr_matrix_times_vector(8, 8, w, v_zero, u, b), std::domain_error);
  std::vector<int> u_zero{0, 3, 5, 7, 7, 9};
  EXPECT_THROW(stan::math::csr_matrix_times_vector(8, 8, w, v_too_long, u_zero, b), std::domain_error);

}
