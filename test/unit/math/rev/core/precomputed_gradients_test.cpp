#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(StanAgradRevInternal, precomputed_gradients) {
  double value;
  std::vector<stan::math::var> vars;
  std::vector<double> gradients;
  stan::math::var x1(2), x2(3);
  stan::math::var y;

  value = 1;
  vars.resize(2);
  vars[0] = x1;
  vars[1] = x2;
  gradients.resize(2);
  gradients[0] = 4;
  gradients[1] = 5;

  EXPECT_NO_THROW(y
                  = stan::math::precomputed_gradients(value, vars, gradients));
  EXPECT_FLOAT_EQ(value, y.val());

  std::vector<double> g;
  EXPECT_NO_THROW(y.grad(vars, g));
  ASSERT_EQ(2U, g.size());
  EXPECT_FLOAT_EQ(gradients[0], g[0]);
  EXPECT_FLOAT_EQ(gradients[1], g[1]);

  stan::math::recover_memory();
}

TEST(StanAgradRevInternal, precomputed_gradients_vari_no_independent_vars) {
  double value = 1;
  std::vector<stan::math::var> vars;
  std::vector<double> gradients;

  stan::math::precomputed_gradients_vari vi(value, vars, gradients);
  EXPECT_FLOAT_EQ(value, vi.val_);
  EXPECT_FLOAT_EQ(0, vi.adj_);
  EXPECT_NO_THROW(vi.chain());
}

TEST(StanAgradRevInternal, precomputed_gradients_vari_mismatched_sizes) {
  double value;
  std::vector<stan::math::var> vars;
  std::vector<double> gradients;

  value = 1;
  vars.resize(1);
  gradients.resize(2);
  EXPECT_THROW(stan::math::precomputed_gradients_vari(value, vars, gradients),
               std::invalid_argument);
}

TEST(StanAgradRevInternal, precomputed_gradients_vari) {
  double value = 1;
  std::vector<stan::math::var> vars;
  stan::math::var x1(2), x2(3);
  vars.push_back(x1);
  vars.push_back(x2);

  std::vector<double> gradients;
  gradients.push_back(4);
  gradients.push_back(5);

  stan::math::precomputed_gradients_vari vi(value, vars, gradients);
  EXPECT_FLOAT_EQ(value, vi.val_);
  EXPECT_FLOAT_EQ(0, vi.adj_);

  EXPECT_NO_THROW(vi.chain())
      << "running vi.chain() with no independent variables";
  EXPECT_FLOAT_EQ(value, vi.val_);
  EXPECT_FLOAT_EQ(0, vi.adj_);
  EXPECT_FLOAT_EQ(0, x1.vi_->adj_);
  EXPECT_FLOAT_EQ(0, x2.vi_->adj_);

  vi.init_dependent();
  EXPECT_NO_THROW(vi.chain())
      << "running vari.chain() with vari initialized as dependent variable";
  EXPECT_FLOAT_EQ(value, vi.val_);
  EXPECT_FLOAT_EQ(1, vi.adj_);
  EXPECT_FLOAT_EQ(gradients[0], x1.vi_->adj_);
  EXPECT_FLOAT_EQ(gradients[1], x2.vi_->adj_);
}

TEST(StanAgradRevInternal, precomputed_gradients_mismatched_sizes) {
  double value;
  std::vector<stan::math::var> vars;
  std::vector<double> gradients;

  value = 1;
  vars.resize(1);
  vars[0] = 0;
  gradients.resize(2);
  gradients[0] = 2;
  gradients[1] = 3;

  EXPECT_THROW(stan::math::precomputed_gradients(value, vars, gradients),
               std::invalid_argument);
  stan::math::recover_memory();
}

TEST(StanAgradRevInternal, precomputed_gradients_containers) {
  double value = 1;
  std::vector<stan::math::var> vars;
  std::vector<double> gradients;
  stan::math::var_value<Eigen::MatrixXd> a(Eigen::MatrixXd::Constant(3, 3, 2));
  std::vector<stan::math::var_value<Eigen::MatrixXd>> b{
      Eigen::MatrixXd::Constant(3, 3, 4), Eigen::MatrixXd::Constant(3, 3, 5)};
  std::tuple<stan::math::var_value<Eigen::MatrixXd>&,
             std::vector<stan::math::var_value<Eigen::MatrixXd>>&>
      ops{a, b};
  Eigen::MatrixXd grad_a = Eigen::MatrixXd::Constant(3, 3, -1);
  std::vector<Eigen::MatrixXd> grad_b{Eigen::MatrixXd::Constant(3, 3, -2),
                                      Eigen::MatrixXd::Constant(3, 3, -3)};
  std::tuple<Eigen::MatrixXd&&, std::vector<Eigen::MatrixXd>&> grads(
      std::move(grad_a), grad_b);

  stan::math::var lp
      = stan::math::precomputed_gradients(value, vars, gradients, ops, grads);
  (2 * lp).grad();
  EXPECT_MATRIX_EQ(a.adj(), Eigen::MatrixXd::Constant(3, 3, -2));
  EXPECT_MATRIX_EQ(b[0].adj(), Eigen::MatrixXd::Constant(3, 3, -4));
  EXPECT_MATRIX_EQ(b[1].adj(), Eigen::MatrixXd::Constant(3, 3, -6));

  stan::math::recover_memory();
}

TEST(StanAgradRevInternal,
     precomputed_gradients_containers_direct_construction) {
  double value = 1;
  std::vector<stan::math::var> vars;
  std::vector<double> gradients;
  stan::math::var_value<Eigen::MatrixXd> a(Eigen::MatrixXd::Constant(3, 3, 2));
  std::vector<stan::math::var_value<Eigen::MatrixXd>> b{
      Eigen::MatrixXd::Constant(3, 3, 4), Eigen::MatrixXd::Constant(3, 3, 5)};
  std::tuple<stan::math::var_value<Eigen::MatrixXd>&,
             std::vector<stan::math::var_value<Eigen::MatrixXd>>&>
      ops{a, b};
  Eigen::MatrixXd grad_a = Eigen::MatrixXd::Constant(3, 3, -1);
  std::vector<Eigen::MatrixXd> grad_b{Eigen::MatrixXd::Constant(3, 3, -2),
                                      Eigen::MatrixXd::Constant(3, 3, -3)};
  std::tuple<Eigen::MatrixXd&&, std::vector<Eigen::MatrixXd>&> grads(
      std::move(grad_a), grad_b);

  stan::math::var lp = new stan::math::precomputed_gradients_vari_template<
      std::tuple<
          stan::arena_t<stan::math::var_value<Eigen::MatrixXd>>,
          stan::arena_t<std::vector<stan::math::var_value<Eigen::MatrixXd>>>>,
      std::tuple<stan::arena_t<Eigen::MatrixXd>,
                 stan::arena_t<std::vector<Eigen::MatrixXd>>>>(
      value, 0, nullptr, nullptr, ops, grads);
  (2 * lp).grad();
  EXPECT_MATRIX_EQ(a.adj(), Eigen::MatrixXd::Constant(3, 3, -2));
  EXPECT_MATRIX_EQ(b[0].adj(), Eigen::MatrixXd::Constant(3, 3, -4));
  EXPECT_MATRIX_EQ(b[1].adj(), Eigen::MatrixXd::Constant(3, 3, -6));

  stan::math::recover_memory();
}

TEST(StanAgradRevInternal, precomputed_gradients_mismatched_containers) {
  double value = 1;
  std::vector<stan::math::var> vars;
  std::vector<double> gradients;
  stan::math::var_value<Eigen::MatrixXd> a(Eigen::MatrixXd::Constant(3, 3, 2));
  std::vector<stan::math::var_value<Eigen::MatrixXd>> b{
      Eigen::MatrixXd::Constant(3, 3, 4), Eigen::MatrixXd::Constant(3, 3, 5)};
  Eigen::MatrixXd grad_a = Eigen::MatrixXd::Constant(3, 2, -1);
  std::vector<Eigen::MatrixXd> grad_b{Eigen::MatrixXd::Constant(3, 2, -2),
                                      Eigen::MatrixXd::Constant(3, 2, -3)};

  EXPECT_THROW(stan::math::precomputed_gradients(value, vars, gradients,
                                                 std::forward_as_tuple(a),
                                                 std::forward_as_tuple(grad_a)),
               std::invalid_argument);
  EXPECT_THROW(stan::math::precomputed_gradients(value, vars, gradients,
                                                 std::forward_as_tuple(b),
                                                 std::forward_as_tuple(grad_b)),
               std::invalid_argument);
  stan::math::recover_memory();
}
