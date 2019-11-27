#include <gtest/gtest.h>
#include <stan/math/prim/arr.hpp>
#include <vector>

TEST(MetaTraits, ScalarSeqViewArray) {
  using stan::scalar_seq_view;
  using std::vector;

  vector<double> v;
  v.push_back(2.2);
  v.push_back(0.0001);
  scalar_seq_view<vector<double> > sv(v);
  EXPECT_FLOAT_EQ(v[0], sv[0]);
  EXPECT_FLOAT_EQ(v[1], sv[1]);

  const vector<double> v_const{2.2, 0.001};
  scalar_seq_view<const vector<double> > sv_const(v_const);
  EXPECT_FLOAT_EQ(v_const[0], sv_const[0]);
  EXPECT_FLOAT_EQ(v_const[1], sv_const[1]);

  const vector<double>& v_const_ref{2.2, 0.001};
  scalar_seq_view<const vector<double> > sv_const_ref(v_const_ref);
  EXPECT_FLOAT_EQ(v_const_ref[0], sv_const_ref[0]);
  EXPECT_FLOAT_EQ(v_const_ref[1], sv_const_ref[1]);

  EXPECT_EQ(v.size(), sv.size());
}
