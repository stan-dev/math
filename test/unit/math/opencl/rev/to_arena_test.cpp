#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, to_arena_matrix_cl_test) {
  Eigen::MatrixXd m(3, 2);
  m << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_cl<double> a(m);

  auto b = stan::math::to_arena(a);
  EXPECT_MATRIX_EQ(stan::math::from_matrix_cl(b),
                   stan::math::from_matrix_cl(a));
  EXPECT_EQ(a.buffer()(), b.buffer()());
  EXPECT_FALSE((std::is_same<decltype(a), decltype(b)>::value));

  auto c = stan::math::to_arena(b);
  EXPECT_TRUE((std::is_same<decltype(b), decltype(c)>::value));
  EXPECT_EQ(b.rows(), c.rows());
  EXPECT_EQ(b.cols(), c.cols());
  EXPECT_EQ(b.buffer()(), c.buffer()());
  stan::math::recover_memory();
}

TEST(AgradRev, to_arena_kg_expression_test) {
  Eigen::MatrixXd m(3, 2);
  m << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_cl<double> a(m);

  auto b = stan::math::to_arena(a + 1);
  EXPECT_MATRIX_EQ(stan::math::from_matrix_cl(b),
                   stan::math::from_matrix_cl(a + 1));
  EXPECT_FALSE((std::is_same<decltype(a), decltype(b)>::value));

  auto c = stan::math::to_arena(b);
  EXPECT_EQ(b.rows(), c.rows());
  EXPECT_EQ(b.cols(), c.cols());
  EXPECT_EQ(b.buffer()(), c.buffer()());
  EXPECT_TRUE((std::is_same<decltype(b), decltype(c)>::value));
  stan::math::recover_memory();
}

TEST(AgradRev, to_arena_var_value_matrix_cl_test) {
  Eigen::MatrixXd val(3, 2);
  val << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd adj(3, 2);
  adj << 4, 5, 6, 7, 8, 9;
  auto* vari = new stan::math::vari_value<stan::math::matrix_cl<double>>(
      stan::math::to_matrix_cl(val));
  vari->adj_ = stan::math::to_matrix_cl(adj);
  stan::math::var_value<stan::math::matrix_cl<double>> a(vari);

  auto b = stan::math::to_arena(a);
  EXPECT_MATRIX_EQ(stan::math::from_matrix_cl(b.val()),
                   stan::math::from_matrix_cl(a.val()));
  EXPECT_MATRIX_EQ(stan::math::from_matrix_cl(b.adj()),
                   stan::math::from_matrix_cl(a.adj()));
  EXPECT_EQ(a.val().buffer()(), b.val().buffer()());
  EXPECT_EQ(a.adj().buffer()(), b.adj().buffer()());
  EXPECT_TRUE((std::is_same<decltype(a), decltype(b)>::value));

  auto c = stan::math::to_arena(b);
  EXPECT_EQ(b.rows(), c.rows());
  EXPECT_EQ(b.cols(), c.cols());
  EXPECT_EQ(b.val().buffer()(), c.val().buffer()());
  EXPECT_EQ(b.adj().buffer()(), c.adj().buffer()());
  EXPECT_TRUE((std::is_same<decltype(b), decltype(c)>::value));
  stan::math::recover_memory();
}

#endif
