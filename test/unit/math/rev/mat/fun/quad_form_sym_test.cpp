#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, quad_form_sym_mat) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::quad_form_sym;

  matrix_v av(4, 4);
  matrix_d ad(4, 4);
  matrix_d bd(4, 2);
  matrix_v bv(4, 2);
  matrix_v res;
  AVEC vars;
  VEC grad;

  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  bv << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;
  av << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  // double-double
  matrix_d resd = quad_form_sym(ad, bd);
  EXPECT_FLOAT_EQ(25433, resd(0, 0));
  EXPECT_FLOAT_EQ(3396, resd(0, 1));
  EXPECT_FLOAT_EQ(3396, resd(1, 0));
  EXPECT_FLOAT_EQ(725, resd(1, 1));
  EXPECT_EQ(resd(1, 0), resd(0, 1));

  // var-double
  res = quad_form_sym(av, bd);
  EXPECT_FLOAT_EQ(25433, res(0, 0).val());
  EXPECT_FLOAT_EQ(3396, res(0, 1).val());
  EXPECT_FLOAT_EQ(3396, res(1, 0).val());
  EXPECT_FLOAT_EQ(725, res(1, 1).val());
  EXPECT_EQ(res(1, 0).val(), res(0, 1).val());

  // double-var
  res = quad_form_sym(ad, bv);
  EXPECT_FLOAT_EQ(25433, res(0, 0).val());
  EXPECT_FLOAT_EQ(3396, res(0, 1).val());
  EXPECT_FLOAT_EQ(3396, res(1, 0).val());
  EXPECT_FLOAT_EQ(725, res(1, 1).val());
  EXPECT_EQ(res(1, 0).val(), res(0, 1).val());

  // var-var
  res = quad_form_sym(av, bv);
  EXPECT_FLOAT_EQ(25433, res(0, 0).val());
  EXPECT_FLOAT_EQ(3396, res(0, 1).val());
  EXPECT_FLOAT_EQ(3396, res(1, 0).val());
  EXPECT_FLOAT_EQ(725, res(1, 1).val());
  EXPECT_EQ(res(1, 0).val(), res(0, 1).val());
}

TEST(AgradRevMatrix, quad_form_sym_mat_grad_vd) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::quad_form_sym;
  using stan::math::sum;

  matrix_v av(4, 4);
  matrix_d ad(4, 4);
  matrix_d bd(4, 2);
  AVAR res;
  AVEC vars;
  VEC grad;
  size_t i, j, pos;

  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;
  av << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  matrix_d dqda(bd * matrix_d::Ones(2, 2) * bd.transpose());

  // var-double
  res = sum(quad_form_sym(av, bd));

  vars.clear();
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      vars.push_back(av(i, j));
  grad = cgradvec(res, vars);
  pos = 0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++, pos++)
      EXPECT_FLOAT_EQ(grad[pos], dqda(i, j));
}

TEST(AgradRevMatrix, quad_form_sym_mat_grad_dv) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::quad_form_sym;
  using stan::math::sum;

  matrix_d ad(4, 4);
  matrix_d bd(4, 2);
  matrix_v bv(4, 2);
  AVAR res;
  AVEC vars;
  VEC grad;
  size_t i, j, pos;

  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  bv << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  matrix_d dqdb((ad * bd + ad.transpose() * bd) * matrix_d::Ones(2, 2));

  // double-var
  res = sum(quad_form_sym(ad, bv));

  vars.clear();
  for (i = 0; i < 4; i++)
    for (j = 0; j < 2; j++)
      vars.push_back(bv(i, j));
  grad = cgradvec(res, vars);
  pos = 0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 2; j++, pos++)
      EXPECT_FLOAT_EQ(grad[pos], dqdb(i, j));
}

TEST(AgradRevMatrix, quad_form_sym_mat_grad_vv) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::quad_form_sym;
  using stan::math::sum;

  matrix_v av(4, 4);
  matrix_d ad(4, 4);
  matrix_d bd(4, 2);
  matrix_v bv(4, 2);
  AVAR res;
  AVEC vars;
  VEC grad;
  size_t i, j, pos;

  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  bv << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;
  av << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  matrix_d dqda(bd * matrix_d::Ones(2, 2) * bd.transpose());
  matrix_d dqdb((ad * bd + ad.transpose() * bd) * matrix_d::Ones(2, 2));

  // var-var
  res = sum(quad_form_sym(av, bv));

  vars.clear();
  for (i = 0; i < 4; i++)
    for (j = 0; j < 2; j++)
      vars.push_back(bv(i, j));
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      vars.push_back(av(i, j));
  grad = cgradvec(res, vars);
  pos = 0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 2; j++, pos++)
      EXPECT_FLOAT_EQ(grad[pos], dqdb(i, j));
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++, pos++)
      EXPECT_FLOAT_EQ(grad[pos], dqda(i, j));
}

TEST(AgradRevMatrix, quad_form_sym_vec) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::quad_form_sym;
  using stan::math::vector_d;
  using stan::math::vector_v;

  matrix_v av(4, 4);
  matrix_d ad(4, 4);
  vector_d bd(4);
  vector_v bv(4);
  stan::math::var res;
  AVEC vars;
  VEC grad;

  bd << 100, 0, -3, 5;
  bv << 100, 0, -3, 5;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;
  av << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  // double-double
  res = quad_form_sym(ad, bd);
  EXPECT_FLOAT_EQ(25433, res.val());

  // var-double
  res = quad_form_sym(av, bd);
  EXPECT_FLOAT_EQ(25433, res.val());

  // double-var
  res = quad_form_sym(ad, bv);
  EXPECT_FLOAT_EQ(25433, res.val());

  // var-var
  res = quad_form_sym(av, bv);
  EXPECT_FLOAT_EQ(25433, res.val());
}

TEST(AgradRevMatrix, quad_form_sym_vec_grad_vd) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::quad_form_sym;
  using stan::math::vector_d;
  using stan::math::vector_v;

  matrix_v av(4, 4);
  matrix_d ad(4, 4);
  vector_d bd(4);
  vector_v bv(4);
  stan::math::var res;
  AVEC vars;
  VEC grad;
  size_t pos, i, j;

  bd << 100, 0, -3, 5;
  bv << 100, 0, -3, 5;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;
  av << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  matrix_d dqda(bd * bd.transpose());
  vector_d dqdb(ad * bd + ad.transpose() * bd);

  // var-double
  res = quad_form_sym(av, bd);

  vars.clear();
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      vars.push_back(av(i, j));
  grad = cgradvec(res, vars);
  for (i = 0, pos = 0; i < 4; i++)
    for (j = 0; j < 4; j++, pos++)
      EXPECT_FLOAT_EQ(grad[pos], dqda(i, j));
}

TEST(AgradRevMatrix, quad_form_sym_vec_grad_dv) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::quad_form_sym;
  using stan::math::vector_d;
  using stan::math::vector_v;

  matrix_v av(4, 4);
  matrix_d ad(4, 4);
  vector_d bd(4);
  vector_v bv(4);
  stan::math::var res;
  AVEC vars;
  VEC grad;
  size_t pos, i;

  bd << 100, 0, -3, 5;
  bv << 100, 0, -3, 5;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;
  av << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  matrix_d dqda(bd * bd.transpose());
  vector_d dqdb(ad * bd + ad.transpose() * bd);

  // double-var
  res = quad_form_sym(ad, bv);

  vars.clear();
  for (i = 0; i < 4; i++)
    vars.push_back(bv[i]);
  grad = cgradvec(res, vars);
  for (i = 0, pos = 0; i < 4; i++, pos++)
    EXPECT_FLOAT_EQ(grad[pos], dqdb[i]);
}

TEST(AgradRevMatrix, quad_form_sym_vec_grad_vv) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::quad_form_sym;
  using stan::math::vector_d;
  using stan::math::vector_v;

  matrix_v av(4, 4);
  matrix_d ad(4, 4);
  vector_d bd(4);
  vector_v bv(4);
  stan::math::var res;
  AVEC vars;
  VEC grad;
  size_t pos, i, j;

  bd << 100, 0, -3, 5;
  bv << 100, 0, -3, 5;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;
  av << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  matrix_d dqda(bd * bd.transpose());
  vector_d dqdb(ad * bd + ad.transpose() * bd);

  // var-var
  res = quad_form_sym(av, bv);

  vars.clear();
  for (size_t i = 0; i < 4; i++)
    vars.push_back(bv[i]);
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      vars.push_back(av(i, j));
  grad = cgradvec(res, vars);
  for (i = 0, pos = 0; i < 4; i++, pos++)
    EXPECT_FLOAT_EQ(grad[pos], dqdb[i]);
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++, pos++)
      EXPECT_FLOAT_EQ(grad[pos], dqda(i, j));
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::to_var;
  stan::math::matrix_d a(4, 4);
  stan::math::matrix_d b1(4, 2);
  stan::math::vector_d b2(4);
  a << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;
  b1 << 100, 10, 0, 1, -3, -3, 5, 2;
  b2 << 100, 10, 0, 1;
  test::check_varis_on_stack(stan::math::quad_form_sym(to_var(a), to_var(b1)));
  test::check_varis_on_stack(stan::math::quad_form_sym(to_var(a), b1));
  test::check_varis_on_stack(stan::math::quad_form_sym(a, to_var(b1)));

  test::check_varis_on_stack(stan::math::quad_form_sym(to_var(a), to_var(b2)));
  test::check_varis_on_stack(stan::math::quad_form_sym(to_var(a), b2));
  test::check_varis_on_stack(stan::math::quad_form_sym(a, to_var(b2)));
}
