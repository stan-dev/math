#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/meta/index_type.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <typeinfo>
#include <vector>

TEST(MathMatrix, Typedefs) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigen_matrix(3, 3);
  Eigen::Matrix<double, Eigen::Dynamic, 1> eigen_vector(3, 1);
  Eigen::Matrix<double, 1, Eigen::Dynamic> eigen_row_vector(1, 3);
  stan::math::index_type<Eigen::Matrix<double,
                                       Eigen::Dynamic,
                                       Eigen::Dynamic> >::type eigen_index;

  stan::math::matrix_d stan_matrix(3, 3);
  stan::math::vector_d stan_vector(3, 1);
  stan::math::row_vector_d stan_row_vector(1, 3);
  stan::math::size_type stan_index;

  std::vector<double> cpp_vector(3);


  EXPECT_EQ(typeid(eigen_index), typeid(stan_index));
  EXPECT_EQ(typeid(eigen_matrix), typeid(stan_matrix));
  EXPECT_EQ(typeid(eigen_vector), typeid(stan_vector));
  EXPECT_EQ(typeid(eigen_row_vector), typeid(stan_row_vector));

  EXPECT_NE(typeid(stan_vector), typeid(cpp_vector));
}
