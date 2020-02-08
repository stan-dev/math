#include <test/unit/math/expect_near_rel.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(testUnitMath, ExpectNearRelScalar) {
  using stan::test::expect_near_rel;
  // zero cases
  expect_near_rel("test A1", 0, 0, 1e-16);
  expect_near_rel("test A2", 0, 0.0, 1e-16);
  expect_near_rel("test A3", 0.0, 0, 1e-16);
  expect_near_rel("test B1", 0, 1e-8, 1e-4);
  expect_near_rel("test B2", 0.0, 1e-8, 1e-4);
  expect_near_rel("test C1", 1e-8, 0, 1e-4);
  expect_near_rel("test C2", 1e-8, 0.0, 1e-4);

  // non-zero examples
  expect_near_rel("test D1", 1, 1, 1e-16);
  expect_near_rel("test D2", 1.0, 1, 1e-16);
  expect_near_rel("test D3", 1, 1.0, 1e-16);
  expect_near_rel("test D4", 1.0, 1.0, 1e-16);
  expect_near_rel("test E1", 1 + 1e-6, 1, 1e-5);
  expect_near_rel("test E2", 1 + 1e-6, 1.0, 1e-5);
  expect_near_rel("test F", 1e4, 1e4 + 1e-6, 1e-8);

  // following fail, but can't unit test through framework
  // expect_near_rel("test G", 0, 1e-6, 1e-8);
  // expect_near_rel("test H", 1e-6, 0, 1e-8);
  // expect_near_rel("test I", 1.0, 1.0 + 1e-6, 1e-8);
  // expect_near_rel("test I", 1.0 + 1e-6, 1.0, 1e-8);
}

TEST(testUnitMath, ExpectNearRelMatrix) {
  using Eigen::Matrix;
  using stan::test::expect_near_rel;
  using stan::test::relative_tolerance;
  typedef Matrix<double, -1, 1> v_t;
  typedef Matrix<double, 1, -1> rv_t;
  typedef Matrix<double, -1, -1> m_t;

  expect_near_rel("test A", v_t(), v_t(), 1e-6);
  expect_near_rel("test A", rv_t(), rv_t(), 1e-6);
  expect_near_rel("test A", m_t(), m_t(), 1e-6);

  v_t b1(2);
  v_t b2(2);
  b1 << 1, 2;
  b2 << 1 + 1e-8, 2 - 1e-8;
  expect_near_rel("test B", b1, b2, 1e-6);

  rv_t c1(2);
  rv_t c2(2);
  c1 << 1, 2;
  c2 << 1 + 1e-8, 2 - 1e-8;
  expect_near_rel("test C", c1, c2, 1e-6);

  m_t d1(2, 3);
  m_t d2(2, 3);
  d1 << 1, 2, 3, 0, 0, 0 - 1e-8;
  d2 << 1 + 1e-8, 2 - 1e-8, 3, 0, 0 + 1e-8, 0;
  expect_near_rel("test D", d1, d2, relative_tolerance(1e-6, 1e-8));

  // these will fail
  // v_t e1(1);
  // v_t e2(2);
  // expect_near_rel("test E", e1, e2, 1e-5);

  // rv_t f1(2);
  // rv_t f2(2);
  // f1 << 1, 2;
  // f2 << 3, 4;
  // expect_near_rel("test F", f1, f2, 1e-2);
}

TEST(testUnitMath, ExpectNearRelVector) {
  using stan::test::expect_near_rel;
  using std::vector;
  typedef std::vector<double> v_t;

  expect_near_rel("test A", v_t{}, v_t{}, 1e-10);

  expect_near_rel("test B", v_t{1, 2, 3}, v_t{1, 2, 3}, 1e-10);

  expect_near_rel("test C", v_t{1, 1, 1}, v_t{1 + 1e-8, 1 - 1e-9, 1}, 1e-6);

  expect_near_rel("test D", v_t{0, 0, 0}, v_t{0, 0 + 9e-9, 0 - 9e-9}, 1e-4);

  // ones after here fail
  // expect_near_rel("test E", v_t{1}, v_t{1, 2}, 1e-6);
  // expect_near_rel("test E", v_t{1, 2}, v_t{}, 1e-6)
}

TEST(testUnitMath, ExpectNearRelVectorNesting) {
  using stan::test::expect_near_rel;
  using stan::test::relative_tolerance;
  using std::vector;
  typedef vector<double> v_t;
  typedef vector<v_t> vv_t;
  typedef vector<vv_t> vvv_t;

  typedef Eigen::Matrix<double, -1, 1> ev_t;
  typedef vector<ev_t> vev_t;

  expect_near_rel("test A", vv_t{}, vv_t{}, 1e-10);

  expect_near_rel("test B", vv_t{v_t{1, 2, 3}, v_t{0, 0, 0}},
                  vv_t{v_t{1, 2, 3}, v_t{0, 0 + 1e-6, 0 - 1e-6}}, 1e-3);

  expect_near_rel(
      "test C",
      vvv_t{vv_t{v_t{1, 2, 3}, v_t{0, 0, 0}}, vv_t{v_t{1, 2, 3}, v_t{0, 0, 0}}},
      vvv_t{vv_t{v_t{1, 2, 3}, v_t{0, 0 + 1e-6, 0 - 1e-6}},
            vv_t{v_t{1, 2, 3}, v_t{0, 0 + 1e-6, 0 - 1e-6}}},
      relative_tolerance(1e-5, 1e-6));

  ev_t d1(3);
  ev_t d2(3);
  d1 << 1, 0, 0;
  d2 << 1 + 1e-8, 0 + 1e-12, 0;
  vev_t e1{d1, d2};
  vev_t e2{d2, d1};
  expect_near_rel("test E", e1, e2, 1e-6);

  // ones after here fail
  // expect_near_rel("test F", v_t{1}, v_t{1, 2}, 1e-6);
  // expect_near_rel("test G", v_t{1, 2}, v_t{}, 1e-6);
  // expect_near_rel("test H", vv_t{v_t{1}, v_t{2}}, vv_t{v_t{1}}, 1e-6);
  // expect_near_rel("test I",
  //                 vev_t{ev_t(3), ev_t(3)},
  //                 vev_t{ev_t(2), ev_t(2)},
  //                 1e-6);
}
