#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <vector>

TEST(AgradMixMatrixEltMultiply, fv_vec_vv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_fv;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);
  fvar<var> c(10.0, 1.0);
  fvar<var> d(100.0, 1.0);

  vector_fv x(2);
  x << a, b;
  vector_fv y(2);
  y << c, d;

  vector_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val());
  EXPECT_FLOAT_EQ(12, z(0).d_.val());
  EXPECT_FLOAT_EQ(105, z(1).d_.val());

  AVEC q = createAVEC(a.val(), b.val(), c.val(), d.val());
  VEC h;
  z(0).val_.grad(q, h);
  EXPECT_FLOAT_EQ(10.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(2.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, fv_vec_vv_2ndDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_fv;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);
  fvar<var> c(10.0, 1.0);
  fvar<var> d(100.0, 1.0);

  vector_fv x(2);
  x << a, b;
  vector_fv y(2);
  y << c, d;

  vector_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val());
  EXPECT_FLOAT_EQ(12, z(0).d_.val());
  EXPECT_FLOAT_EQ(105, z(1).d_.val());

  AVEC q = createAVEC(a.val(), b.val(), c.val(), d.val());
  VEC h;
  z(0).d_.grad(q, h);
  EXPECT_FLOAT_EQ(1.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(1.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, fv_vec_vd_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_fv;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);

  vector_fv x(2);
  x << a, b;
  vector_d y(2);
  y << 10, 100;

  vector_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val());
  EXPECT_FLOAT_EQ(10, z(0).d_.val());
  EXPECT_FLOAT_EQ(100, z(1).d_.val());

  AVEC q = createAVEC(a.val(), b.val());
  VEC h;
  z(0).val_.grad(q, h);
  EXPECT_FLOAT_EQ(10.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, fv_vec_vd_2ndDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_fv;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);

  vector_fv x(2);
  x << a, b;
  vector_d y(2);
  y << 10, 100;

  vector_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val());
  EXPECT_FLOAT_EQ(10, z(0).d_.val());
  EXPECT_FLOAT_EQ(100, z(1).d_.val());

  AVEC q = createAVEC(a.val(), b.val());
  VEC h;
  z(0).d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, fv_vec_dv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_fv;

  fvar<var> c(10.0, 1.0);
  fvar<var> d(100.0, 1.0);

  vector_d x(2);
  x << 2, 5;
  vector_fv y(2);
  y << c, d;

  vector_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val());
  EXPECT_FLOAT_EQ(2, z(0).d_.val());
  EXPECT_FLOAT_EQ(5, z(1).d_.val());

  AVEC q = createAVEC(c.val(), d.val());
  VEC h;
  z(0).val_.grad(q, h);
  EXPECT_FLOAT_EQ(2.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, fv_vec_dv_2ndDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_fv;

  fvar<var> c(10.0, 1.0);
  fvar<var> d(100.0, 1.0);

  vector_d x(2);
  x << 2, 5;
  vector_fv y(2);
  y << c, d;

  vector_fv z = elt_multiply(x, y);

  AVEC q = createAVEC(c.val(), d.val());
  VEC h;
  z(0).d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, fv_row_vec_vv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_fv;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);
  fvar<var> c(10.0, 1.0);
  fvar<var> d(100.0, 1.0);

  row_vector_fv x(2);
  x << a, b;
  row_vector_fv y(2);
  y << c, d;

  row_vector_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val());
  EXPECT_FLOAT_EQ(12, z(0).d_.val());
  EXPECT_FLOAT_EQ(105, z(1).d_.val());

  AVEC q = createAVEC(a.val(), b.val(), c.val(), d.val());
  VEC h;
  z(0).val_.grad(q, h);
  EXPECT_FLOAT_EQ(10.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(2.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, fv_row_vec_vv_2ndDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_fv;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);
  fvar<var> c(10.0, 1.0);
  fvar<var> d(100.0, 1.0);

  row_vector_fv x(2);
  x << a, b;
  row_vector_fv y(2);
  y << c, d;

  row_vector_fv z = elt_multiply(x, y);

  AVEC q = createAVEC(a.val(), b.val(), c.val(), d.val());
  VEC h;
  z(0).d_.grad(q, h);
  EXPECT_FLOAT_EQ(1.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(1.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, fv_row_vec_vd_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_fv;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);

  row_vector_fv x(2);
  x << a, b;
  row_vector_d y(2);
  y << 10, 100;

  row_vector_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val());
  EXPECT_FLOAT_EQ(10, z(0).d_.val());
  EXPECT_FLOAT_EQ(100, z(1).d_.val());

  AVEC q = createAVEC(a.val(), b.val());
  VEC h;
  z(0).val_.grad(q, h);
  EXPECT_FLOAT_EQ(10.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, fv_row_vec_vd_2ndDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_fv;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);

  row_vector_fv x(2);
  x << a, b;
  row_vector_d y(2);
  y << 10, 100;

  row_vector_fv z = elt_multiply(x, y);

  AVEC q = createAVEC(a.val(), b.val());
  VEC h;
  z(0).d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, fv_row_vec_dv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_fv;
  using stan::math::var;

  fvar<var> c(10.0, 1.0);
  fvar<var> d(100.0, 1.0);

  row_vector_d x(2);
  x << 2, 5;
  row_vector_fv y(2);
  y << c, d;

  row_vector_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val());
  EXPECT_FLOAT_EQ(2, z(0).d_.val());
  EXPECT_FLOAT_EQ(5, z(1).d_.val());

  AVEC q = createAVEC(c.val(), d.val());
  VEC h;
  z(0).val_.grad(q, h);
  EXPECT_FLOAT_EQ(2.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, fv_row_vec_dv_2ndDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_fv;
  using stan::math::var;

  fvar<var> c(10.0, 1.0);
  fvar<var> d(100.0, 1.0);

  row_vector_d x(2);
  x << 2, 5;
  row_vector_fv y(2);
  y << c, d;

  row_vector_fv z = elt_multiply(x, y);

  AVEC q = createAVEC(c.val(), d.val());
  VEC h;
  z(0).d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}

TEST(AgradMixMatrixEltMultiply, fv_matrix_vv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_fv;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);
  fvar<var> c(6.0, 1.0);
  fvar<var> d(9.0, 1.0);
  fvar<var> e(13.0, 1.0);
  fvar<var> f(29.0, 1.0);
  fvar<var> g(10.0, 1.0);
  fvar<var> h(100.0, 1.0);
  fvar<var> i(1000.0, 1.0);
  fvar<var> j(10000.0, 1.0);
  fvar<var> k(100000.0, 1.0);
  fvar<var> l(1000000.0, 1.0);

  matrix_fv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_fv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val());
  EXPECT_FLOAT_EQ(12, z(0, 0).d_.val());
  EXPECT_FLOAT_EQ(105, z(0, 1).d_.val());
  EXPECT_FLOAT_EQ(1000029, z(1, 2).d_.val());

  std::vector<var> q;
  q.push_back(a.val());
  q.push_back(b.val());
  q.push_back(c.val());
  q.push_back(d.val());
  q.push_back(e.val());
  q.push_back(f.val());
  q.push_back(g.val());
  q.push_back(h.val());
  q.push_back(i.val());
  q.push_back(j.val());
  q.push_back(k.val());
  q.push_back(l.val());
  VEC hh;
  z(0).val_.grad(q, hh);
  EXPECT_FLOAT_EQ(10, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
  EXPECT_FLOAT_EQ(2, hh[6]);
  EXPECT_FLOAT_EQ(0.0, hh[7]);
  EXPECT_FLOAT_EQ(0.0, hh[8]);
  EXPECT_FLOAT_EQ(0.0, hh[9]);
  EXPECT_FLOAT_EQ(0.0, hh[10]);
  EXPECT_FLOAT_EQ(0.0, hh[11]);
}
TEST(AgradMixMatrixEltMultiply, fv_matrix_vv_2ndDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_fv;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);
  fvar<var> c(6.0, 1.0);
  fvar<var> d(9.0, 1.0);
  fvar<var> e(13.0, 1.0);
  fvar<var> f(29.0, 1.0);
  fvar<var> g(10.0, 1.0);
  fvar<var> h(100.0, 1.0);
  fvar<var> i(1000.0, 1.0);
  fvar<var> j(10000.0, 1.0);
  fvar<var> k(100000.0, 1.0);
  fvar<var> l(1000000.0, 1.0);

  matrix_fv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_fv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val());
  EXPECT_FLOAT_EQ(12, z(0, 0).d_.val());
  EXPECT_FLOAT_EQ(105, z(0, 1).d_.val());
  EXPECT_FLOAT_EQ(1000029, z(1, 2).d_.val());

  std::vector<var> q;
  q.push_back(a.val());
  q.push_back(b.val());
  q.push_back(c.val());
  q.push_back(d.val());
  q.push_back(e.val());
  q.push_back(f.val());
  q.push_back(g.val());
  q.push_back(h.val());
  q.push_back(i.val());
  q.push_back(j.val());
  q.push_back(k.val());
  q.push_back(l.val());
  VEC hh;
  z(0).d_.grad(q, hh);
  EXPECT_FLOAT_EQ(1, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
  EXPECT_FLOAT_EQ(1, hh[6]);
  EXPECT_FLOAT_EQ(0.0, hh[7]);
  EXPECT_FLOAT_EQ(0.0, hh[8]);
  EXPECT_FLOAT_EQ(0.0, hh[9]);
  EXPECT_FLOAT_EQ(0.0, hh[10]);
  EXPECT_FLOAT_EQ(0.0, hh[11]);
}
TEST(AgradMixMatrixEltMultiply, fv_matrix_vd_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_fv;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);
  fvar<var> c(6.0, 1.0);
  fvar<var> d(9.0, 1.0);
  fvar<var> e(13.0, 1.0);
  fvar<var> f(29.0, 1.0);

  matrix_fv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_d y(2, 3);
  y << 10, 100, 1000, 10000, 100000, 1000000;

  matrix_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val());
  EXPECT_FLOAT_EQ(10, z(0, 0).d_.val());
  EXPECT_FLOAT_EQ(100, z(0, 1).d_.val());
  EXPECT_FLOAT_EQ(1000000, z(1, 2).d_.val());

  AVEC q = createAVEC(a.val(), b.val(), c.val(), d.val(), e.val(), f.val());
  VEC h;
  z(0).val_.grad(q, h);
  EXPECT_FLOAT_EQ(10, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
  EXPECT_FLOAT_EQ(0.0, h[4]);
  EXPECT_FLOAT_EQ(0.0, h[5]);
}
TEST(AgradMixMatrixEltMultiply, fv_matrix_vd_2ndDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_fv;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(5.0, 1.0);
  fvar<var> c(6.0, 1.0);
  fvar<var> d(9.0, 1.0);
  fvar<var> e(13.0, 1.0);
  fvar<var> f(29.0, 1.0);

  matrix_fv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_d y(2, 3);
  y << 10, 100, 1000, 10000, 100000, 1000000;

  matrix_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val());
  EXPECT_FLOAT_EQ(10, z(0, 0).d_.val());
  EXPECT_FLOAT_EQ(100, z(0, 1).d_.val());
  EXPECT_FLOAT_EQ(1000000, z(1, 2).d_.val());

  AVEC q = createAVEC(a.val(), b.val(), c.val(), d.val(), e.val(), f.val());
  VEC h;
  z(0).d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
  EXPECT_FLOAT_EQ(0.0, h[4]);
  EXPECT_FLOAT_EQ(0.0, h[5]);
}
TEST(AgradMixMatrixEltMultiply, fv_matrix_dv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_fv;
  using stan::math::var;

  fvar<var> g(10.0, 1.0);
  fvar<var> h(100.0, 1.0);
  fvar<var> i(1000.0, 1.0);
  fvar<var> j(10000.0, 1.0);
  fvar<var> k(100000.0, 1.0);
  fvar<var> l(1000000.0, 1.0);

  matrix_d x(2, 3);
  x << 2, 5, 6, 9, 13, 29;
  matrix_fv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val());
  EXPECT_FLOAT_EQ(2, z(0, 0).d_.val());
  EXPECT_FLOAT_EQ(5, z(0, 1).d_.val());
  EXPECT_FLOAT_EQ(29, z(1, 2).d_.val());

  AVEC q = createAVEC(g.val(), h.val(), i.val(), j.val(), k.val(), l.val());
  VEC hh;
  z(0).val_.grad(q, hh);
  EXPECT_FLOAT_EQ(2, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
}
TEST(AgradMixMatrixEltMultiply, fv_matrix_dv_2ndDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_fv;
  using stan::math::var;

  fvar<var> g(10.0, 1.0);
  fvar<var> h(100.0, 1.0);
  fvar<var> i(1000.0, 1.0);
  fvar<var> j(10000.0, 1.0);
  fvar<var> k(100000.0, 1.0);
  fvar<var> l(1000000.0, 1.0);

  matrix_d x(2, 3);
  x << 2, 5, 6, 9, 13, 29;
  matrix_fv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_fv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val());
  EXPECT_FLOAT_EQ(2, z(0, 0).d_.val());
  EXPECT_FLOAT_EQ(5, z(0, 1).d_.val());
  EXPECT_FLOAT_EQ(29, z(1, 2).d_.val());

  AVEC q = createAVEC(g.val(), h.val(), i.val(), j.val(), k.val(), l.val());
  VEC hh;
  z(0).d_.grad(q, hh);
  EXPECT_FLOAT_EQ(0, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_vv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_ffv;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  vector_ffv x(2);
  x << a, b;
  vector_ffv y(2);
  y << c, d;

  vector_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val().val());
  EXPECT_FLOAT_EQ(12, z(0).d_.val().val());
  EXPECT_FLOAT_EQ(105, z(1).d_.val().val());

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  z(0).val_.val().grad(q, h);
  EXPECT_FLOAT_EQ(10.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(2.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_vv_2ndDeriv_1) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_ffv;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  vector_ffv x(2);
  x << a, b;
  vector_ffv y(2);
  y << c, d;

  vector_ffv z = elt_multiply(x, y);

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  z(0).val().d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_vv_2ndDeriv_2) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_ffv;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  vector_ffv x(2);
  x << a, b;
  vector_ffv y(2);
  y << c, d;

  vector_ffv z = elt_multiply(x, y);

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  z(0).d_.val().grad(q, h);
  EXPECT_FLOAT_EQ(1.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(1.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_vv_3rdDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_ffv;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  vector_ffv x(2);
  x << a, b;
  vector_ffv y(2);
  y << c, d;

  vector_ffv z = elt_multiply(x, y);

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  z(0).d_.d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_vd_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);

  vector_ffv x(2);
  x << a, b;
  vector_d y(2);
  y << 10, 100;

  vector_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val().val());
  EXPECT_FLOAT_EQ(10, z(0).d_.val().val());
  EXPECT_FLOAT_EQ(100, z(1).d_.val().val());

  AVEC q = createAVEC(a.val().val(), b.val().val());
  VEC h;
  z(0).val_.val().grad(q, h);
  EXPECT_FLOAT_EQ(10.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_vd_2ndDeriv_1) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);

  vector_ffv x(2);
  x << a, b;
  vector_d y(2);
  y << 10, 100;

  vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(a.val().val(), b.val().val());
  VEC h;
  z(0).val().d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_vd_2ndDeriv_2) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);

  vector_ffv x(2);
  x << a, b;
  vector_d y(2);
  y << 10, 100;

  vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(a.val().val(), b.val().val());
  VEC h;
  z(0).d_.val().grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_vd_3rdDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  a.val_.d_ = 1.0;
  b.val_.d_ = 1.0;

  vector_ffv x(2);
  x << a, b;
  vector_d y(2);
  y << 10, 100;

  vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(a.val().val(), b.val().val());
  VEC h;
  z(0).d_.d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_dv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;

  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  vector_d x(2);
  x << 2, 5;
  vector_ffv y(2);
  y << c, d;

  vector_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val().val());
  EXPECT_FLOAT_EQ(2, z(0).d_.val().val());
  EXPECT_FLOAT_EQ(5, z(1).d_.val().val());

  AVEC q = createAVEC(c.val().val(), d.val().val());
  VEC h;
  z(0).val_.val().grad(q, h);
  EXPECT_FLOAT_EQ(2.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_dv_2ndDeriv_1) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;

  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  vector_d x(2);
  x << 2, 5;
  vector_ffv y(2);
  y << c, d;

  vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(c.val().val(), d.val().val());
  VEC h;
  z(0).val().d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_dv_2ndDeriv_2) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;

  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  vector_d x(2);
  x << 2, 5;
  vector_ffv y(2);
  y << c, d;

  vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(c.val().val(), d.val().val());
  VEC h;
  z(0).d_.val().grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_vec_dv_3rdDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;

  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);
  c.val_.d_ = 1.0;
  d.val_.d_ = 1.0;

  vector_d x(2);
  x << 2, 5;
  vector_ffv y(2);
  y << c, d;

  vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(c.val().val(), d.val().val());
  VEC h;
  z(0).d_.d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_vv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  row_vector_ffv x(2);
  x << a, b;
  row_vector_ffv y(2);
  y << c, d;

  row_vector_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val().val());
  EXPECT_FLOAT_EQ(12, z(0).d_.val().val());
  EXPECT_FLOAT_EQ(105, z(1).d_.val().val());

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  z(0).val_.val().grad(q, h);
  EXPECT_FLOAT_EQ(10.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(2.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_vv_2ndDeriv_1) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  row_vector_ffv x(2);
  x << a, b;
  row_vector_ffv y(2);
  y << c, d;

  row_vector_ffv z = elt_multiply(x, y);

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  z(0).val().d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_vv_2ndDeriv_2) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  row_vector_ffv x(2);
  x << a, b;
  row_vector_ffv y(2);
  y << c, d;

  row_vector_ffv z = elt_multiply(x, y);

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  z(0).d_.val().grad(q, h);
  EXPECT_FLOAT_EQ(1.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(1.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_vv_3rdDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  row_vector_ffv x(2);
  x << a, b;
  row_vector_ffv y(2);
  y << c, d;

  row_vector_ffv z = elt_multiply(x, y);

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  z(0).d_.d_.grad(q, h);
  EXPECT_FLOAT_EQ(0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_vd_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);

  row_vector_ffv x(2);
  x << a, b;
  row_vector_d y(2);
  y << 10, 100;

  row_vector_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val().val());
  EXPECT_FLOAT_EQ(10, z(0).d_.val().val());
  EXPECT_FLOAT_EQ(100, z(1).d_.val().val());

  AVEC q = createAVEC(a.val().val(), b.val().val());
  VEC h;
  z(0).val_.val().grad(q, h);
  EXPECT_FLOAT_EQ(10.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_vd_2ndDeriv_1) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);

  row_vector_ffv x(2);
  x << a, b;
  row_vector_d y(2);
  y << 10, 100;

  row_vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(a.val().val(), b.val().val());
  VEC h;
  z(0).val().d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_vd_2ndDeriv_2) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);

  row_vector_ffv x(2);
  x << a, b;
  row_vector_d y(2);
  y << 10, 100;

  row_vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(a.val().val(), b.val().val());
  VEC h;
  z(0).d_.val().grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_vd_3rdDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  a.val_.d_ = 1.0;
  b.val_.d_ = 1.0;

  row_vector_ffv x(2);
  x << a, b;
  row_vector_d y(2);
  y << 10, 100;

  row_vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(a.val().val(), b.val().val());
  VEC h;
  z(0).d_.d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_dv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  row_vector_d x(2);
  x << 2, 5;
  row_vector_ffv y(2);
  y << c, d;

  row_vector_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(1).val_.val().val());
  EXPECT_FLOAT_EQ(2, z(0).d_.val().val());
  EXPECT_FLOAT_EQ(5, z(1).d_.val().val());

  AVEC q = createAVEC(c.val().val(), d.val().val());
  VEC h;
  z(0).val_.val().grad(q, h);
  EXPECT_FLOAT_EQ(2.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_dv_2ndDeriv_1) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  row_vector_d x(2);
  x << 2, 5;
  row_vector_ffv y(2);
  y << c, d;

  row_vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(c.val().val(), d.val().val());
  VEC h;
  z(0).val().d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_dv_2ndDeriv_2) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);

  row_vector_d x(2);
  x << 2, 5;
  row_vector_ffv y(2);
  y << c, d;

  row_vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(c.val().val(), d.val().val());
  VEC h;
  z(0).d_.val().grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_row_vec_dv_3rdDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::row_vector_d;
  using stan::math::row_vector_ffv;
  using stan::math::var;

  fvar<fvar<var> > c(10.0, 1.0);
  fvar<fvar<var> > d(100.0, 1.0);
  c.val_.d_ = 1.0;
  d.val_.d_ = 1.0;

  row_vector_d x(2);
  x << 2, 5;
  row_vector_ffv y(2);
  y << c, d;

  row_vector_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(c.val().val(), d.val().val());
  VEC h;
  z(0).d_.d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_vv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(6.0, 1.0);
  fvar<fvar<var> > d(9.0, 1.0);
  fvar<fvar<var> > e(13.0, 1.0);
  fvar<fvar<var> > f(29.0, 1.0);
  fvar<fvar<var> > g(10.0, 1.0);
  fvar<fvar<var> > h(100.0, 1.0);
  fvar<fvar<var> > i(1000.0, 1.0);
  fvar<fvar<var> > j(10000.0, 1.0);
  fvar<fvar<var> > k(100000.0, 1.0);
  fvar<fvar<var> > l(1000000.0, 1.0);

  matrix_ffv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_ffv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val().val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val().val());
  EXPECT_FLOAT_EQ(12, z(0, 0).d_.val().val());
  EXPECT_FLOAT_EQ(105, z(0, 1).d_.val().val());
  EXPECT_FLOAT_EQ(1000029, z(1, 2).d_.val().val());

  std::vector<var> q;
  q.push_back(a.val().val());
  q.push_back(b.val().val());
  q.push_back(c.val().val());
  q.push_back(d.val().val());
  q.push_back(e.val().val());
  q.push_back(f.val().val());
  q.push_back(g.val().val());
  q.push_back(h.val().val());
  q.push_back(i.val().val());
  q.push_back(j.val().val());
  q.push_back(k.val().val());
  q.push_back(l.val().val());
  VEC hh;
  z(0).val_.val().grad(q, hh);
  EXPECT_FLOAT_EQ(10, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
  EXPECT_FLOAT_EQ(2, hh[6]);
  EXPECT_FLOAT_EQ(0.0, hh[7]);
  EXPECT_FLOAT_EQ(0.0, hh[8]);
  EXPECT_FLOAT_EQ(0.0, hh[9]);
  EXPECT_FLOAT_EQ(0.0, hh[10]);
  EXPECT_FLOAT_EQ(0.0, hh[11]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_vv_2ndDeriv_1) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(6.0, 1.0);
  fvar<fvar<var> > d(9.0, 1.0);
  fvar<fvar<var> > e(13.0, 1.0);
  fvar<fvar<var> > f(29.0, 1.0);
  fvar<fvar<var> > g(10.0, 1.0);
  fvar<fvar<var> > h(100.0, 1.0);
  fvar<fvar<var> > i(1000.0, 1.0);
  fvar<fvar<var> > j(10000.0, 1.0);
  fvar<fvar<var> > k(100000.0, 1.0);
  fvar<fvar<var> > l(1000000.0, 1.0);

  matrix_ffv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_ffv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_ffv z = elt_multiply(x, y);

  std::vector<var> q;
  q.push_back(a.val().val());
  q.push_back(b.val().val());
  q.push_back(c.val().val());
  q.push_back(d.val().val());
  q.push_back(e.val().val());
  q.push_back(f.val().val());
  q.push_back(g.val().val());
  q.push_back(h.val().val());
  q.push_back(i.val().val());
  q.push_back(j.val().val());
  q.push_back(k.val().val());
  q.push_back(l.val().val());
  VEC hh;
  z(0).val().d_.grad(q, hh);
  EXPECT_FLOAT_EQ(0, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
  EXPECT_FLOAT_EQ(0, hh[6]);
  EXPECT_FLOAT_EQ(0.0, hh[7]);
  EXPECT_FLOAT_EQ(0.0, hh[8]);
  EXPECT_FLOAT_EQ(0.0, hh[9]);
  EXPECT_FLOAT_EQ(0.0, hh[10]);
  EXPECT_FLOAT_EQ(0.0, hh[11]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_vv_2ndDeriv_2) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(6.0, 1.0);
  fvar<fvar<var> > d(9.0, 1.0);
  fvar<fvar<var> > e(13.0, 1.0);
  fvar<fvar<var> > f(29.0, 1.0);
  fvar<fvar<var> > g(10.0, 1.0);
  fvar<fvar<var> > h(100.0, 1.0);
  fvar<fvar<var> > i(1000.0, 1.0);
  fvar<fvar<var> > j(10000.0, 1.0);
  fvar<fvar<var> > k(100000.0, 1.0);
  fvar<fvar<var> > l(1000000.0, 1.0);

  matrix_ffv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_ffv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_ffv z = elt_multiply(x, y);

  std::vector<var> q;
  q.push_back(a.val().val());
  q.push_back(b.val().val());
  q.push_back(c.val().val());
  q.push_back(d.val().val());
  q.push_back(e.val().val());
  q.push_back(f.val().val());
  q.push_back(g.val().val());
  q.push_back(h.val().val());
  q.push_back(i.val().val());
  q.push_back(j.val().val());
  q.push_back(k.val().val());
  q.push_back(l.val().val());
  VEC hh;
  z(0).d_.val().grad(q, hh);
  EXPECT_FLOAT_EQ(1, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
  EXPECT_FLOAT_EQ(1, hh[6]);
  EXPECT_FLOAT_EQ(0.0, hh[7]);
  EXPECT_FLOAT_EQ(0.0, hh[8]);
  EXPECT_FLOAT_EQ(0.0, hh[9]);
  EXPECT_FLOAT_EQ(0.0, hh[10]);
  EXPECT_FLOAT_EQ(0.0, hh[11]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_vv_3rdDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(6.0, 1.0);
  fvar<fvar<var> > d(9.0, 1.0);
  fvar<fvar<var> > e(13.0, 1.0);
  fvar<fvar<var> > f(29.0, 1.0);
  fvar<fvar<var> > g(10.0, 1.0);
  fvar<fvar<var> > h(100.0, 1.0);
  fvar<fvar<var> > i(1000.0, 1.0);
  fvar<fvar<var> > j(10000.0, 1.0);
  fvar<fvar<var> > k(100000.0, 1.0);
  fvar<fvar<var> > l(1000000.0, 1.0);

  matrix_ffv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_ffv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_ffv z = elt_multiply(x, y);

  std::vector<var> q;
  q.push_back(a.val().val());
  q.push_back(b.val().val());
  q.push_back(c.val().val());
  q.push_back(d.val().val());
  q.push_back(e.val().val());
  q.push_back(f.val().val());
  q.push_back(g.val().val());
  q.push_back(h.val().val());
  q.push_back(i.val().val());
  q.push_back(j.val().val());
  q.push_back(k.val().val());
  q.push_back(l.val().val());
  VEC hh;
  z(0).d_.d_.grad(q, hh);
  EXPECT_FLOAT_EQ(0, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
  EXPECT_FLOAT_EQ(0, hh[6]);
  EXPECT_FLOAT_EQ(0.0, hh[7]);
  EXPECT_FLOAT_EQ(0.0, hh[8]);
  EXPECT_FLOAT_EQ(0.0, hh[9]);
  EXPECT_FLOAT_EQ(0.0, hh[10]);
  EXPECT_FLOAT_EQ(0.0, hh[11]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_vd_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(6.0, 1.0);
  fvar<fvar<var> > d(9.0, 1.0);
  fvar<fvar<var> > e(13.0, 1.0);
  fvar<fvar<var> > f(29.0, 1.0);

  matrix_ffv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_d y(2, 3);
  y << 10, 100, 1000, 10000, 100000, 1000000;

  matrix_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val().val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val().val());
  EXPECT_FLOAT_EQ(10, z(0, 0).d_.val().val());
  EXPECT_FLOAT_EQ(100, z(0, 1).d_.val().val());
  EXPECT_FLOAT_EQ(1000000, z(1, 2).d_.val().val());

  AVEC q = createAVEC(a.val().val(), b.val().val(), c.val().val(),
                      d.val().val(), e.val().val(), f.val().val());
  VEC h;
  z(0).val_.val().grad(q, h);
  EXPECT_FLOAT_EQ(10, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
  EXPECT_FLOAT_EQ(0.0, h[4]);
  EXPECT_FLOAT_EQ(0.0, h[5]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_vd_2ndDeriv_1) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(6.0, 1.0);
  fvar<fvar<var> > d(9.0, 1.0);
  fvar<fvar<var> > e(13.0, 1.0);
  fvar<fvar<var> > f(29.0, 1.0);

  matrix_ffv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_d y(2, 3);
  y << 10, 100, 1000, 10000, 100000, 1000000;

  matrix_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val().val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val().val());
  EXPECT_FLOAT_EQ(10, z(0, 0).d_.val().val());
  EXPECT_FLOAT_EQ(100, z(0, 1).d_.val().val());
  EXPECT_FLOAT_EQ(1000000, z(1, 2).d_.val().val());

  AVEC q = createAVEC(a.val().val(), b.val().val(), c.val().val(),
                      d.val().val(), e.val().val(), f.val().val());
  VEC h;
  z(0).val().d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
  EXPECT_FLOAT_EQ(0.0, h[4]);
  EXPECT_FLOAT_EQ(0.0, h[5]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_vd_2ndDeriv_2) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(6.0, 1.0);
  fvar<fvar<var> > d(9.0, 1.0);
  fvar<fvar<var> > e(13.0, 1.0);
  fvar<fvar<var> > f(29.0, 1.0);

  matrix_ffv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_d y(2, 3);
  y << 10, 100, 1000, 10000, 100000, 1000000;

  matrix_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val().val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val().val());
  EXPECT_FLOAT_EQ(10, z(0, 0).d_.val().val());
  EXPECT_FLOAT_EQ(100, z(0, 1).d_.val().val());
  EXPECT_FLOAT_EQ(1000000, z(1, 2).d_.val().val());

  AVEC q = createAVEC(a.val().val(), b.val().val(), c.val().val(),
                      d.val().val(), e.val().val(), f.val().val());
  VEC h;
  z(0).d_.val().grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
  EXPECT_FLOAT_EQ(0.0, h[4]);
  EXPECT_FLOAT_EQ(0.0, h[5]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_vd_3rdDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(2.0, 1.0);
  fvar<fvar<var> > b(5.0, 1.0);
  fvar<fvar<var> > c(6.0, 1.0);
  fvar<fvar<var> > d(9.0, 1.0);
  fvar<fvar<var> > e(13.0, 1.0);
  fvar<fvar<var> > f(29.0, 1.0);
  a.val_.d_ = 1.0;
  b.val_.d_ = 1.0;
  c.val_.d_ = 1.0;
  d.val_.d_ = 1.0;
  e.val_.d_ = 1.0;
  f.val_.d_ = 1.0;

  matrix_ffv x(2, 3);
  x << a, b, c, d, e, f;
  matrix_d y(2, 3);
  y << 10, 100, 1000, 10000, 100000, 1000000;

  matrix_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val().val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val().val());
  EXPECT_FLOAT_EQ(10, z(0, 0).d_.val().val());
  EXPECT_FLOAT_EQ(100, z(0, 1).d_.val().val());
  EXPECT_FLOAT_EQ(1000000, z(1, 2).d_.val().val());

  AVEC q = createAVEC(a.val().val(), b.val().val(), c.val().val(),
                      d.val().val(), e.val().val(), f.val().val());
  VEC h;
  z(0).d_.d_.grad(q, h);
  EXPECT_FLOAT_EQ(0.0, h[0]);
  EXPECT_FLOAT_EQ(0.0, h[1]);
  EXPECT_FLOAT_EQ(0.0, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
  EXPECT_FLOAT_EQ(0.0, h[4]);
  EXPECT_FLOAT_EQ(0.0, h[5]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_dv_1stDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > g(10.0, 1.0);
  fvar<fvar<var> > h(100.0, 1.0);
  fvar<fvar<var> > i(1000.0, 1.0);
  fvar<fvar<var> > j(10000.0, 1.0);
  fvar<fvar<var> > k(100000.0, 1.0);
  fvar<fvar<var> > l(1000000.0, 1.0);

  matrix_d x(2, 3);
  x << 2, 5, 6, 9, 13, 29;
  matrix_ffv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_ffv z = elt_multiply(x, y);
  EXPECT_FLOAT_EQ(20.0, z(0, 0).val_.val().val());
  EXPECT_FLOAT_EQ(500.0, z(0, 1).val_.val().val());
  EXPECT_FLOAT_EQ(29000000.0, z(1, 2).val_.val().val());
  EXPECT_FLOAT_EQ(2, z(0, 0).d_.val().val());
  EXPECT_FLOAT_EQ(5, z(0, 1).d_.val().val());
  EXPECT_FLOAT_EQ(29, z(1, 2).d_.val().val());

  AVEC q = createAVEC(g.val().val(), h.val().val(), i.val().val(),
                      j.val().val(), k.val().val(), l.val().val());
  VEC hh;
  z(0).val_.val().grad(q, hh);
  EXPECT_FLOAT_EQ(2, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_dv_2ndDeriv_1) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > g(10.0, 1.0);
  fvar<fvar<var> > h(100.0, 1.0);
  fvar<fvar<var> > i(1000.0, 1.0);
  fvar<fvar<var> > j(10000.0, 1.0);
  fvar<fvar<var> > k(100000.0, 1.0);
  fvar<fvar<var> > l(1000000.0, 1.0);

  matrix_d x(2, 3);
  x << 2, 5, 6, 9, 13, 29;
  matrix_ffv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(g.val().val(), h.val().val(), i.val().val(),
                      j.val().val(), k.val().val(), l.val().val());
  VEC hh;
  z(0).val().d_.grad(q, hh);
  EXPECT_FLOAT_EQ(0, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_dv_2ndDeriv_2) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > g(10.0, 1.0);
  fvar<fvar<var> > h(100.0, 1.0);
  fvar<fvar<var> > i(1000.0, 1.0);
  fvar<fvar<var> > j(10000.0, 1.0);
  fvar<fvar<var> > k(100000.0, 1.0);
  fvar<fvar<var> > l(1000000.0, 1.0);

  matrix_d x(2, 3);
  x << 2, 5, 6, 9, 13, 29;
  matrix_ffv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(g.val().val(), h.val().val(), i.val().val(),
                      j.val().val(), k.val().val(), l.val().val());
  VEC hh;
  z(0).d_.val().grad(q, hh);
  EXPECT_FLOAT_EQ(0, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
}
TEST(AgradMixMatrixEltMultiply, ffv_matrix_dv_3rdDeriv) {
  using stan::math::elt_multiply;
  using stan::math::fvar;
  using stan::math::matrix_d;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > g(10.0, 1.0);
  fvar<fvar<var> > h(100.0, 1.0);
  fvar<fvar<var> > i(1000.0, 1.0);
  fvar<fvar<var> > j(10000.0, 1.0);
  fvar<fvar<var> > k(100000.0, 1.0);
  fvar<fvar<var> > l(1000000.0, 1.0);
  g.val_.d_ = 1.0;
  h.val_.d_ = 1.0;
  i.val_.d_ = 1.0;
  j.val_.d_ = 1.0;
  k.val_.d_ = 1.0;
  l.val_.d_ = 1.0;

  matrix_d x(2, 3);
  x << 2, 5, 6, 9, 13, 29;
  matrix_ffv y(2, 3);
  y << g, h, i, j, k, l;

  matrix_ffv z = elt_multiply(x, y);

  AVEC q = createAVEC(g.val().val(), h.val().val(), i.val().val(),
                      j.val().val(), k.val().val(), l.val().val());
  VEC hh;
  z(0).d_.d_.grad(q, hh);
  EXPECT_FLOAT_EQ(0, hh[0]);
  EXPECT_FLOAT_EQ(0.0, hh[1]);
  EXPECT_FLOAT_EQ(0.0, hh[2]);
  EXPECT_FLOAT_EQ(0.0, hh[3]);
  EXPECT_FLOAT_EQ(0.0, hh[4]);
  EXPECT_FLOAT_EQ(0.0, hh[5]);
}
