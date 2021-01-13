#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaRev, ScalarSeqViewvar) {
  using stan::scalar_seq_view;
  using stan::math::var;
  var d = 10;
  scalar_seq_view<var> sv(d);
  EXPECT_FLOAT_EQ(d.val(), sv.val(0));
  EXPECT_FLOAT_EQ(d.val(), sv.val(12));

  const var d_const = 10;
  scalar_seq_view<const var> sv_const(d_const);
  EXPECT_FLOAT_EQ(d_const.val(), sv_const.val(0));
  EXPECT_FLOAT_EQ(d_const.val(), sv_const.val(12));

  const var& d_const_ref = 10;
  scalar_seq_view<const var&> sv_const_ref(d_const);
  EXPECT_FLOAT_EQ(d_const_ref.val(), sv_const_ref.val(0));
  EXPECT_FLOAT_EQ(d_const_ref.val(), sv_const_ref.val(12));

  EXPECT_EQ(1, sv.size());
}

TEST(MathMetaRev, ScalarSeqViewArrayVarVal) {
  using stan::scalar_seq_view;
  using stan::math::var;
  using std::vector;
  vector<var> v;
  v.push_back(2.2);
  v.push_back(0.0001);
  scalar_seq_view<vector<var>> sv(v);
  EXPECT_FLOAT_EQ(v[0].val(), sv.val(0));
  EXPECT_FLOAT_EQ(v[1].val(), sv.val(1));

  const vector<var> v_const{2.2, 0.001};
  scalar_seq_view<const vector<var>> sv_const(v_const);
  EXPECT_FLOAT_EQ(v_const[0].val(), sv_const.val(0));
  EXPECT_FLOAT_EQ(v_const[1].val(), sv_const.val(1));

  const vector<var>& v_const_ref{2.2, 0.001};
  scalar_seq_view<const vector<var>> sv_const_ref(v_const_ref);
  EXPECT_FLOAT_EQ(v_const_ref[0].val(), sv_const_ref.val(0));
  EXPECT_FLOAT_EQ(v_const_ref[1].val(), sv_const_ref.val(1));

  EXPECT_EQ(v.size(), sv.size());
}

template <typename C>
void expect_scalar_seq_view_value(const C& v) {
  using stan::scalar_seq_view;
  scalar_seq_view<C> sv(v);
  EXPECT_FLOAT_EQ(v.val()(0), sv.val(0));
  EXPECT_FLOAT_EQ(v.val()(1), sv.val(1));
  EXPECT_FLOAT_EQ(v.val()(2), sv.val(2));
  EXPECT_FLOAT_EQ(v.val()(3), sv.val(3));

  EXPECT_EQ(v.size(), sv.size());
}

template <typename C>
void expect_scalar_seq_view_adjoints(const C& v) {
  using stan::scalar_seq_view;
  scalar_seq_view<C> sv(v);
  std::vector<stan::math::var> stdv(sv.size());
  for (size_t i = 0; i < sv.size(); ++i) {
    stdv[i] = sv[i];
  }

  for (size_t i = 0; i < sv.size(); ++i) {
    stan::math::set_zero_all_adjoints();
    stdv[i].grad();
    EXPECT_EQ(1.0, v.adj()(i));
    for (size_t j = 0; j < sv.size(); ++j) {
      if (j != i) {
        EXPECT_EQ(0.0, v.adj()(j));
      }
    }
  }
}

TEST(MathMetaRev, ScalarSeqViewVectorVar) {
  using stan::math::var;
  Eigen::Matrix<var, -1, 1> A = Eigen::VectorXd(4);
  expect_scalar_seq_view_value(A);
  expect_scalar_seq_view_adjoints(A);
}

TEST(MathMetaRev, ScalarSeqViewRowVectorVar) {
  using stan::math::var;
  Eigen::Matrix<var, 1, -1> A = Eigen::RowVectorXd(4);
  expect_scalar_seq_view_value(A);
  expect_scalar_seq_view_adjoints(A);
}

TEST(MathMetaRev, VarScalarSeqViewVector) {
  using stan::math::var_value;
  var_value<Eigen::Matrix<double, -1, 1>> A = Eigen::VectorXd(4);
  expect_scalar_seq_view_value(A);
  expect_scalar_seq_view_adjoints(A);
}

TEST(MathMetaRev, VarScalarSeqViewRowVector) {
  using stan::math::var_value;
  var_value<Eigen::Matrix<double, 1, -1>> A = Eigen::RowVectorXd(4);
  expect_scalar_seq_view_value(A);
  expect_scalar_seq_view_adjoints(A);
}
