#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRev, rev_matrix_type_test) {
  using v_matrix = stan::math::var_value<Eigen::MatrixXd>;
  using v_vector = stan::math::var_value<Eigen::VectorXd>;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::vector_d;
  using stan::math::vector_v;
  using t1 = stan::rev_matrix_t<-1, -1, v_vector>;
  EXPECT_SAME_TYPE(t1, v_matrix);
  using t2 = stan::rev_matrix_t<-1, -1, vector_v>;
  EXPECT_SAME_TYPE(t2, matrix_v);
  using t3 = stan::rev_matrix_t<-1, 1, vector_d>;
  EXPECT_SAME_TYPE(t3, vector_d);

  using t4 = stan::rev_matrix_t<-1, 1, vector_d, vector_v, matrix_d>;
  EXPECT_SAME_TYPE(t4, vector_v);
  using t5 = stan::rev_matrix_t<-1, -1, vector_v, v_vector, matrix_d>;
  EXPECT_SAME_TYPE(t5, v_matrix);

  using t6 = stan::rev_matrix_t<-1, 1, vector_d, vector_d, matrix_d>;
  EXPECT_SAME_TYPE(t6, vector_d);
  using t7 = stan::rev_matrix_t<-1, -1, vector_v, vector_d, matrix_d>;
  EXPECT_SAME_TYPE(t7, matrix_v);
}
