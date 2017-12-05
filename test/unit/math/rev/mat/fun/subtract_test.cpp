#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, subtract_scalar) {
  using stan::math::subtract;
  using stan::math::matrix_v;

  matrix_v v(2, 2);
  v << 1, 2, 3, 4;
  matrix_v result;

  result = subtract(2.0, v);
  EXPECT_FLOAT_EQ(1.0, result(0, 0).val());
  EXPECT_FLOAT_EQ(0.0, result(0, 1).val());
  EXPECT_FLOAT_EQ(-1.0, result(1, 0).val());
  EXPECT_FLOAT_EQ(-2.0, result(1, 1).val());

  result = subtract(v, 2.0);
  EXPECT_FLOAT_EQ(-1.0, result(0, 0).val());
  EXPECT_FLOAT_EQ(0.0, result(0, 1).val());
  EXPECT_FLOAT_EQ(1.0, result(1, 0).val());
  EXPECT_FLOAT_EQ(2.0, result(1, 1).val());
}

TEST(AgradRevMatrix, subtract_vector_vector) {
  using stan::math::subtract;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d expected_output(5);
  vector_v output;
  vector_d output_d;
  vector_d vd_1(5), vd_2(5);
  vector_v vv_1(5), vv_2(5);

  vd_1 << 0, 2, -6, 10, 6;
  vv_1 << 0, 2, -6, 10, 6;
  vd_2 << 2, 3, 4, 5, 6;
  vv_2 << 2, 3, 4, 5, 6;

  expected_output << -2, -1, -10, 5, 0;

  output_d = subtract(vd_1, vd_2);
  EXPECT_FLOAT_EQ(expected_output(0), output_d(0));
  EXPECT_FLOAT_EQ(expected_output(1), output_d(1));
  EXPECT_FLOAT_EQ(expected_output(2), output_d(2));
  EXPECT_FLOAT_EQ(expected_output(3), output_d(3));
  EXPECT_FLOAT_EQ(expected_output(4), output_d(4));

  output = subtract(vv_1, vd_2);
  EXPECT_FLOAT_EQ(expected_output(0), output(0).val());
  EXPECT_FLOAT_EQ(expected_output(1), output(1).val());
  EXPECT_FLOAT_EQ(expected_output(2), output(2).val());
  EXPECT_FLOAT_EQ(expected_output(3), output(3).val());
  EXPECT_FLOAT_EQ(expected_output(4), output(4).val());

  output = subtract(vd_1, vv_2);
  EXPECT_FLOAT_EQ(expected_output(0), output(0).val());
  EXPECT_FLOAT_EQ(expected_output(1), output(1).val());
  EXPECT_FLOAT_EQ(expected_output(2), output(2).val());
  EXPECT_FLOAT_EQ(expected_output(3), output(3).val());
  EXPECT_FLOAT_EQ(expected_output(4), output(4).val());

  output = subtract(vv_1, vv_2);
  EXPECT_FLOAT_EQ(expected_output(0), output(0).val());
  EXPECT_FLOAT_EQ(expected_output(1), output(1).val());
  EXPECT_FLOAT_EQ(expected_output(2), output(2).val());
  EXPECT_FLOAT_EQ(expected_output(3), output(3).val());
  EXPECT_FLOAT_EQ(expected_output(4), output(4).val());
}
TEST(AgradRevMatrix, subtract_vector_vector_exception) {
  using stan::math::subtract;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d1(5), d2(1);
  vector_v v1(5), v2(1);

  vector_v output;
  EXPECT_THROW(subtract(d1, d2), std::invalid_argument);
  EXPECT_THROW(subtract(v1, d2), std::invalid_argument);
  EXPECT_THROW(subtract(d1, v2), std::invalid_argument);
  EXPECT_THROW(subtract(v1, v2), std::invalid_argument);
}
TEST(AgradRevMatrix, subtract_rowvector_rowvector) {
  using stan::math::subtract;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;

  row_vector_d expected_output(5);
  row_vector_d  output_d;
  row_vector_v  output;
  row_vector_d rvd_1(5), rvd_2(5);
  row_vector_v rvv_1(5), rvv_2(5);

  rvd_1 << 0, 2, -6, 10, 6;
  rvv_1 << 0, 2, -6, 10, 6;
  rvd_2 << 2, 3, 4, 5, 6;
  rvv_2 << 2, 3, 4, 5, 6;

  expected_output << -2, -1, -10, 5, 0;

  output_d = subtract(rvd_1, rvd_2);
  EXPECT_FLOAT_EQ(expected_output(0), output_d(0));
  EXPECT_FLOAT_EQ(expected_output(1), output_d(1));
  EXPECT_FLOAT_EQ(expected_output(2), output_d(2));
  EXPECT_FLOAT_EQ(expected_output(3), output_d(3));
  EXPECT_FLOAT_EQ(expected_output(4), output_d(4));

  output = subtract(rvv_1, rvd_2);
  EXPECT_FLOAT_EQ(expected_output(0), output(0).val());
  EXPECT_FLOAT_EQ(expected_output(1), output(1).val());
  EXPECT_FLOAT_EQ(expected_output(2), output(2).val());
  EXPECT_FLOAT_EQ(expected_output(3), output(3).val());
  EXPECT_FLOAT_EQ(expected_output(4), output(4).val());

  output = subtract(rvd_1, rvv_2);
  EXPECT_FLOAT_EQ(expected_output(0), output(0).val());
  EXPECT_FLOAT_EQ(expected_output(1), output(1).val());
  EXPECT_FLOAT_EQ(expected_output(2), output(2).val());
  EXPECT_FLOAT_EQ(expected_output(3), output(3).val());
  EXPECT_FLOAT_EQ(expected_output(4), output(4).val());

  output = subtract(rvv_1, rvv_2);
  EXPECT_FLOAT_EQ(expected_output(0), output(0).val());
  EXPECT_FLOAT_EQ(expected_output(1), output(1).val());
  EXPECT_FLOAT_EQ(expected_output(2), output(2).val());
  EXPECT_FLOAT_EQ(expected_output(3), output(3).val());
  EXPECT_FLOAT_EQ(expected_output(4), output(4).val());
}
TEST(AgradRevMatrix, subtract_rowvector_rowvector_exception) {
  using stan::math::subtract;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;

  row_vector_d d1(5), d2(2);
  row_vector_v v1(5), v2(2);

  row_vector_v output;
  EXPECT_THROW(subtract(d1, d2), std::invalid_argument);
  EXPECT_THROW(subtract(d1, v2), std::invalid_argument);
  EXPECT_THROW(subtract(v1, d2), std::invalid_argument);
  EXPECT_THROW(subtract(v1, v2), std::invalid_argument);
}
TEST(AgradRevMatrix, subtract_matrix_matrix) {
  using stan::math::subtract;
  using stan::math::matrix_d;
  using stan::math::matrix_v;

  matrix_d expected_output(2, 2);
  matrix_v output;
  matrix_d md_1(2, 2), md_2(2, 2);
  matrix_v mv_1(2, 2), mv_2(2, 2);
  matrix_d md_mis(2, 3);
  matrix_v mv_mis(1, 1);

  md_1 << -10, 1, 10, 0;
  mv_1 << -10, 1, 10, 0;
  md_2 << 10, -10, 1, 2;
  mv_2 << 10, -10, 1, 2;

  expected_output << -20, 11, 9, -2;

  matrix_d output_d = subtract(md_1, md_2);
  EXPECT_FLOAT_EQ(expected_output(0, 0), output_d(0, 0));
  EXPECT_FLOAT_EQ(expected_output(0, 1), output_d(0, 1));
  EXPECT_FLOAT_EQ(expected_output(1, 0), output_d(1, 0));
  EXPECT_FLOAT_EQ(expected_output(1, 1), output_d(1, 1));

  output = subtract(mv_1, md_2);
  EXPECT_FLOAT_EQ(expected_output(0, 0), output(0, 0).val());
  EXPECT_FLOAT_EQ(expected_output(0, 1), output(0, 1).val());
  EXPECT_FLOAT_EQ(expected_output(1, 0), output(1, 0).val());
  EXPECT_FLOAT_EQ(expected_output(1, 1), output(1, 1).val());

  output = subtract(md_1, mv_2);
  EXPECT_FLOAT_EQ(expected_output(0, 0), output(0, 0).val());
  EXPECT_FLOAT_EQ(expected_output(0, 1), output(0, 1).val());
  EXPECT_FLOAT_EQ(expected_output(1, 0), output(1, 0).val());
  EXPECT_FLOAT_EQ(expected_output(1, 1), output(1, 1).val());

  output = subtract(mv_1, mv_2);
  EXPECT_FLOAT_EQ(expected_output(0, 0), output(0, 0).val());
  EXPECT_FLOAT_EQ(expected_output(0, 1), output(0, 1).val());
  EXPECT_FLOAT_EQ(expected_output(1, 0), output(1, 0).val());
  EXPECT_FLOAT_EQ(expected_output(1, 1), output(1, 1).val());
}
TEST(AgradRevMatrix, subtract_matrix_matrix_exception) {
  using stan::math::subtract;
  using stan::math::matrix_d;
  using stan::math::matrix_v;

  matrix_d d1(2, 2), d2(1, 2);
  matrix_v v1(2, 2), v2(1, 2);

  EXPECT_THROW(subtract(d1, d2), std::invalid_argument);
  EXPECT_THROW(subtract(d1, v2), std::invalid_argument);
  EXPECT_THROW(subtract(v1, d2), std::invalid_argument);
  EXPECT_THROW(subtract(v1, v2), std::invalid_argument);
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::value_of;
  stan::math::matrix_v m(2, 2);
  m << 1, 2, 3, 4;

  stan::math::vector_v v(3);
  v << 1, 2, 3;

  stan::math::row_vector_v rv(2);
  rv << 1, 2;

  stan::math::var s(2);

  test::check_varis_on_stack(stan::math::subtract(m, m));
  test::check_varis_on_stack(stan::math::subtract(m, value_of(m)));
  test::check_varis_on_stack(stan::math::subtract(value_of(m), m));

  test::check_varis_on_stack(stan::math::subtract(v, v));
  test::check_varis_on_stack(stan::math::subtract(v, value_of(v)));
  test::check_varis_on_stack(stan::math::subtract(value_of(v), v));

  test::check_varis_on_stack(stan::math::subtract(rv, rv));
  test::check_varis_on_stack(stan::math::subtract(rv, value_of(rv)));
  test::check_varis_on_stack(stan::math::subtract(value_of(rv), rv));

  test::check_varis_on_stack(stan::math::subtract(m, s));
  test::check_varis_on_stack(stan::math::subtract(m, value_of(s)));
  test::check_varis_on_stack(stan::math::subtract(value_of(m), s));

  test::check_varis_on_stack(stan::math::subtract(v, s));
  test::check_varis_on_stack(stan::math::subtract(v, value_of(s)));
  test::check_varis_on_stack(stan::math::subtract(value_of(v), s));

  test::check_varis_on_stack(stan::math::subtract(rv, s));
  test::check_varis_on_stack(stan::math::subtract(rv, value_of(s)));
  test::check_varis_on_stack(stan::math::subtract(value_of(rv), s));

  test::check_varis_on_stack(stan::math::subtract(s, m));
  test::check_varis_on_stack(stan::math::subtract(s, value_of(m)));
  test::check_varis_on_stack(stan::math::subtract(value_of(s), m));

  test::check_varis_on_stack(stan::math::subtract(s, v));
  test::check_varis_on_stack(stan::math::subtract(s, value_of(v)));
  test::check_varis_on_stack(stan::math::subtract(value_of(s), v));

  test::check_varis_on_stack(stan::math::subtract(s, rv));
  test::check_varis_on_stack(stan::math::subtract(s, value_of(rv)));
  test::check_varis_on_stack(stan::math::subtract(value_of(s), rv));
}
