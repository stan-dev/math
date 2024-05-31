#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, csr_matrix_times_vector_vals) {
  using stan::math::csr_matrix_times_vector;
  std::vector<int> v{1, 2, 3, 1, 2};
  std::vector<int> u{1, 2, 3, 4, 5, 6};
  Eigen::VectorXd b(5);
  b << 1, 2, 3, 4, 5;
  Eigen::VectorXd w(5);
  w << -0.67082, 0.5, -0.223607, -0.223607, -0.5;
  auto dbl_res = csr_matrix_times_vector(5, 5, w, v, u, b);
  using stan::math::var;
  Eigen::Matrix<var, -1, 1> b_var(b);
  Eigen::Matrix<var, -1, 1> w_var(w);
  auto var_res_w = csr_matrix_times_vector(5, 5, w_var, v, u, b);
  EXPECT_MATRIX_EQ(dbl_res, var_res_w);
  auto var_res_b = csr_matrix_times_vector(5, 5, w, v, u, b_var);
  EXPECT_MATRIX_EQ(dbl_res, var_res_b);
  auto var_res_bw = csr_matrix_times_vector(5, 5, w_var, v, u, b_var);
  EXPECT_MATRIX_EQ(dbl_res, var_res_bw);
}

TEST(MathMixMatFun, csr_matrix_times_vector1) {
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
}

TEST(MathMixMatFun, csr_matrix_times_vector2) {
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

  stan::test::expect_ad_matvar(f, w, b);
}
