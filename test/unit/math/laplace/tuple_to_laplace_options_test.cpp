#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>

TEST(laplace_utils, tuple_to_laplace_options) {
  using stan::math::laplace_options_user_supplied;
  using stan::math::internal::tuple_to_laplace_options;

  auto ops = std::make_tuple(Eigen::VectorXd::Zero(3), 1e-6, 100, 1, 2, 0);
  auto laplace_opts = tuple_to_laplace_options(ops);
  EXPECT_EQ(laplace_opts.hessian_block_size, 1);
  EXPECT_EQ(laplace_opts.solver, 1);
  EXPECT_EQ(laplace_opts.tolerance, 1e-6);
  EXPECT_EQ(laplace_opts.max_num_steps, 100);
  EXPECT_EQ(laplace_opts.line_search.max_iterations, 2);
  EXPECT_EQ(laplace_opts.allow_fallthrough, false);
  EXPECT_EQ(laplace_opts.theta_0, Eigen::VectorXd::Zero(3));
  static_assert(std::is_same_v<decltype(laplace_opts), laplace_options_user_supplied>);
}

TEST(laplace_utils, tuple_to_laplace_options_move) {
  using stan::math::laplace_options_user_supplied;
  using stan::math::internal::tuple_to_laplace_options;

  auto ops = std::make_tuple(Eigen::VectorXd::Zero(3), 1e-6, 100, 1, 2, 1);
  auto laplace_opts = tuple_to_laplace_options(std::move(ops));
  EXPECT_EQ(laplace_opts.hessian_block_size, 1);
  EXPECT_EQ(laplace_opts.solver, 1);
  EXPECT_EQ(laplace_opts.tolerance, 1e-6);
  EXPECT_EQ(laplace_opts.max_num_steps, 100);
  EXPECT_EQ(laplace_opts.line_search.max_iterations, 2);
  EXPECT_EQ(laplace_opts.allow_fallthrough, true);
  EXPECT_EQ(laplace_opts.theta_0, Eigen::VectorXd::Zero(3));
  static_assert(std::is_same_v<decltype(laplace_opts), laplace_options_user_supplied>);
}
