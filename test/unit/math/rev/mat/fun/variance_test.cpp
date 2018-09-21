#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <vector>

TEST(AgradRevMatrix, varianceZeroBoundaryCase) {
  using stan::math::var;
  using stan::math::variance;
  using std::vector;

  vector<var> y(3, 1.7);
  var f = variance(y);
  EXPECT_FLOAT_EQ(0.0, f.val());

  vector<double> g;
  f.grad(y, g);
  EXPECT_EQ(y.size(), g.size());
  for (size_t i = 0; i < g.size(); ++i)
    EXPECT_FLOAT_EQ(0.0, g[i]);
}

TEST(AgradRevMatrix, variance_vector) {
  using stan::math::variance;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d(1);
  d << 12.9;
  EXPECT_FLOAT_EQ(0.0, variance(d));

  vector_d d1(6);
  vector_v v1(6);

  d1 << 1, 2, 3, 4, 5, 6;
  v1 << 1, 2, 3, 4, 5, 6;

  EXPECT_FLOAT_EQ(17.5 / 5.0, variance(d1));

  EXPECT_FLOAT_EQ(17.5 / 5.0, variance(v1).val());

  d1.resize(1);
  v1.resize(1);
  EXPECT_FLOAT_EQ(0.0, variance(d1));
  EXPECT_FLOAT_EQ(0.0, variance(v1).val());
}

TEST(AgradRevMatrix, variance_avoid_precision_loss) {
  using stan::math::variance;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_v v(4);
  v << 1.0, 2.0, 4.0, 5.0;
  stan::math::var output = variance(v);
  output.grad();

  EXPECT_FLOAT_EQ(10.0 / 3.0, output.val());

  for (int i = 0; i < v.size(); ++i) {
    EXPECT_FLOAT_EQ(2.0 * (v(i).val() - 3.0) / 3.0, v(i).adj());
  }
}

TEST(AgradRevMatrix, variance_vector_exception) {
  using stan::math::variance;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d1;
  vector_v v1;
  EXPECT_THROW(variance(d1), std::invalid_argument);
  EXPECT_THROW(variance(v1), std::invalid_argument);
}
TEST(AgradRevMatrix, variance_rowvector) {
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;
  using stan::math::variance;

  row_vector_d d(1);
  d << 12.9;
  EXPECT_FLOAT_EQ(0.0, variance(d));

  row_vector_d d1(6);
  row_vector_v v1(6);

  d1 << 1, 2, 3, 4, 5, 6;
  v1 << 1, 2, 3, 4, 5, 6;

  EXPECT_FLOAT_EQ(17.5 / 5.0, variance(d1));

  EXPECT_FLOAT_EQ(17.5 / 5.0, variance(v1).val());

  d1.resize(1);
  v1.resize(1);
  EXPECT_FLOAT_EQ(0.0, variance(d1));
  EXPECT_FLOAT_EQ(0.0, variance(v1).val());
}
TEST(AgradRevMatrix, variance_rowvector_exception) {
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;
  using stan::math::variance;

  row_vector_d d1;
  row_vector_v v1;
  EXPECT_THROW(variance(d1), std::invalid_argument);
  EXPECT_THROW(variance(v1), std::invalid_argument);
}
TEST(AgradRevMatrix, variance_matrix) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::variance;

  matrix_d m(1, 1);
  m << 12.9;
  EXPECT_FLOAT_EQ(0.0, variance(m));

  matrix_d d1(2, 3);
  matrix_v v1(2, 3);

  d1 << 1, 2, 3, 4, 5, 6;
  v1 << 1, 2, 3, 4, 5, 6;

  EXPECT_FLOAT_EQ(17.5 / 5.0, variance(d1));

  EXPECT_FLOAT_EQ(17.5 / 5.0, variance(v1).val());

  d1.resize(1, 1);
  v1.resize(1, 1);
  EXPECT_FLOAT_EQ(0.0, variance(d1));
  EXPECT_FLOAT_EQ(0.0, variance(v1).val());
}
TEST(AgradRevMatrix, variance_matrix_exception) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::variance;

  matrix_d d1;
  matrix_v v1;
  EXPECT_THROW(variance(d1), std::invalid_argument);
  EXPECT_THROW(variance(v1), std::invalid_argument);

  d1.resize(0, 1);
  v1.resize(0, 1);
  EXPECT_THROW(variance(d1), std::invalid_argument);
  EXPECT_THROW(variance(v1), std::invalid_argument);

  d1.resize(1, 0);
  v1.resize(1, 0);
  EXPECT_THROW(variance(d1), std::invalid_argument);
  EXPECT_THROW(variance(v1), std::invalid_argument);
}
TEST(AgradRevMatrix, varianceStdVector) {
  // should use arg-dep lookup
  using stan::math::variance;

  AVEC y1 = createAVEC(0.5, 2.0, 3.5);
  AVAR f1 = variance(y1);
  VEC grad1 = cgrad(f1, y1[0], y1[1], y1[2]);
  // save before cleaned out
  double f1_val = f1.val();

  AVEC y2 = createAVEC(0.5, 2.0, 3.5);
  AVAR mean2 = (y2[0] + y2[1] + y2[2]) / 3.0;
  AVAR sum_sq_diff_2 = (y2[0] - mean2) * (y2[0] - mean2)
                       + (y2[1] - mean2) * (y2[1] - mean2)
                       + (y2[2] - mean2) * (y2[2] - mean2);
  AVAR f2 = sum_sq_diff_2 / (3 - 1);

  EXPECT_EQ(f2.val(), f1_val);

  VEC grad2 = cgrad(f2, y2[0], y2[1], y2[2]);

  EXPECT_EQ(3U, grad1.size());
  EXPECT_EQ(3U, grad2.size());
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(grad2[i], grad1[i]);
  }
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::vector_v v1(6);
  v1 << 1, 2, 3, 4, 5, 6;
  stan::math::row_vector_v v2(6);
  v2 << 1, 2, 3, 4, 5, 6;

  test::check_varis_on_stack(stan::math::variance(v1));
  test::check_varis_on_stack(stan::math::variance(v2));
}
