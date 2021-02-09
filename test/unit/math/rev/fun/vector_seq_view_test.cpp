#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <vector>

template <typename T>
void test_vec_seq_var(const T& m1) {
  using stan::vector_seq_view;
  using stan::math::var;
  vector_seq_view<T> vsv(m1);
  EXPECT_FLOAT_EQ(m1.val()(1), vsv.val(12)(1));
  EXPECT_EQ(vsv.size(), 1);

  std::vector<stan::plain_type_t<T>> v;
  v.push_back(m1);
  v.push_back(m1.reverse());
  vector_seq_view<std::vector<stan::plain_type_t<T>>> vsv_vec(v);
  EXPECT_FLOAT_EQ(m1.val()(0), vsv_vec.val(0)[0]);
  EXPECT_FLOAT_EQ(m1.val()(0), vsv_vec.val(1)[3]);
  EXPECT_FLOAT_EQ(m1.val()(2), vsv_vec.val(1)[1]);
  EXPECT_FLOAT_EQ(m1.val()(2), vsv_vec.val(0)[2]);
  EXPECT_EQ(vsv_vec.size(), 2);
}
TEST(MathMetaRev, VectorSeqViewVar) {
  Eigen::Matrix<double, -1, 1> values = Eigen::Matrix<double, -1, 1>::Random(4);
  Eigen::Matrix<stan::math::var, -1, 1> A = values;
  test_vec_seq_var(A);
  test_vec_seq_var(A.transpose());
  stan::math::var_value<Eigen::Matrix<double, -1, 1>> B = values;
  test_vec_seq_var(B);
  test_vec_seq_var(B.transpose());
}
