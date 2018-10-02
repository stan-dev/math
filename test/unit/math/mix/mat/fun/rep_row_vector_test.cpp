#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>

using stan::math::var;
TEST(AgradMixMatrixRepRowVector, fv_rowvector) {
  using stan::math::fvar;
  using stan::math::rep_row_vector;
  using stan::math::row_vector_fv;
  fvar<var> a;
  a.val_ = 3.0;
  a.d_ = 2.0;
  row_vector_fv output;
  output = rep_row_vector(a, 4);

  EXPECT_EQ(3, output(0).val_.val());
  EXPECT_EQ(3, output(1).val_.val());
  EXPECT_EQ(3, output(2).val_.val());
  EXPECT_EQ(3, output(3).val_.val());
  EXPECT_EQ(2, output(0).d_.val());
  EXPECT_EQ(2, output(1).d_.val());
  EXPECT_EQ(2, output(2).d_.val());
  EXPECT_EQ(2, output(3).d_.val());
}

TEST(AgradMixMatrixRepRowVector, fv_rowvector_exception) {
  using stan::math::fvar;
  using stan::math::rep_row_vector;
  using stan::math::row_vector_fv;
  fvar<var> a;
  a.val_ = 3.0;
  a.d_ = 2.0;

  EXPECT_THROW(rep_row_vector(a, -2), std::domain_error);
}
TEST(AgradMixMatrixRepRowVector, ffv_rowvector) {
  using stan::math::fvar;
  using stan::math::rep_row_vector;
  using stan::math::row_vector_ffv;
  fvar<fvar<var> > a;
  a.val_ = 3.0;
  a.d_ = 2.0;
  row_vector_ffv output;
  output = rep_row_vector(a, 4);

  EXPECT_EQ(3, output(0).val_.val().val());
  EXPECT_EQ(3, output(1).val_.val().val());
  EXPECT_EQ(3, output(2).val_.val().val());
  EXPECT_EQ(3, output(3).val_.val().val());
  EXPECT_EQ(2, output(0).d_.val().val());
  EXPECT_EQ(2, output(1).d_.val().val());
  EXPECT_EQ(2, output(2).d_.val().val());
  EXPECT_EQ(2, output(3).d_.val().val());
}

TEST(AgradMixMatrixRepRowVector, ffv_rowvector_exception) {
  using stan::math::fvar;
  using stan::math::rep_row_vector;
  using stan::math::row_vector_ffv;
  fvar<fvar<var> > a;
  a.val_ = 3.0;
  a.d_ = 2.0;

  EXPECT_THROW(rep_row_vector(a, -2), std::domain_error);
}
