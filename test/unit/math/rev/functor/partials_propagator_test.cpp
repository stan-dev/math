#include <stan/math/rev/functor/partials_propagator.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradPartialsVari, PartialsPropagatorScal) {
  using stan::math::make_partials_propagator;
  using stan::math::var;

  double d1;
  auto o3 = stan::math::make_partials_propagator(d1);
  EXPECT_EQ(2, sizeof(o3));

  var v1 = var(0.0);

  std::vector<var> v_stdvec;
  v_stdvec.push_back(v1);

  auto o4 = stan::math::make_partials_propagator(v1);
  stan::math::edge<0>(o4).partials_[0] += 10.0;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
}

TEST(AgradPartialsVari, PartialsPropagatorVec) {
  using stan::math::make_partials_propagator;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d_vec(4);
  auto o3 = stan::math::make_partials_propagator(d_vec);
  EXPECT_EQ(2, sizeof(o3));

  vector_v v_vec(4);
  var v1 = var(0.0);
  var v2 = var(1.0);
  var v3 = var(2.0);
  var v4 = var(3.0);
  v_vec << v1, v2, v3, v4;

  std::vector<var> v_stdvec;
  v_stdvec.push_back(v1);
  v_stdvec.push_back(v2);
  v_stdvec.push_back(v3);
  v_stdvec.push_back(v4);

  auto o4 = stan::math::make_partials_propagator(v_vec);
  stan::math::edge<0>(o4).partials_[0] += 10.0;
  stan::math::edge<0>(o4).partials_[1] += 20.0;
  stan::math::edge<0>(o4).partials_[2] += 30.0;
  stan::math::edge<0>(o4).partials_[3] += 40.0;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
  EXPECT_FLOAT_EQ(20.0, grad[1]);
  EXPECT_FLOAT_EQ(30.0, grad[2]);
  EXPECT_FLOAT_EQ(40.0, grad[3]);
}

TEST(AgradPartialsVari, PartialsPropagatorStdVec) {
  using stan::math::make_partials_propagator;
  using stan::math::var;

  std::vector<double> d_vec(4);
  auto o3 = stan::math::make_partials_propagator(d_vec);
  EXPECT_EQ(2, sizeof(o3));

  std::vector<var> v_vec;
  var v1 = var(0.0);
  var v2 = var(1.0);
  var v3 = var(2.0);
  var v4 = var(3.0);
  v_vec.push_back(v1);
  v_vec.push_back(v2);
  v_vec.push_back(v3);
  v_vec.push_back(v4);

  auto o4 = stan::math::make_partials_propagator(v_vec);
  stan::math::edge<0>(o4).partials_[0] += 10.0;
  stan::math::edge<0>(o4).partials_[1] += 20.0;
  stan::math::edge<0>(o4).partials_[2] += 30.0;
  stan::math::edge<0>(o4).partials_[3] += 40.0;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_vec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
  EXPECT_FLOAT_EQ(20.0, grad[1]);
  EXPECT_FLOAT_EQ(30.0, grad[2]);
  EXPECT_FLOAT_EQ(40.0, grad[3]);
}

TEST(AgradPartialsVari, PartialsPropagatorMat) {
  using stan::math::make_partials_propagator;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::var;

  matrix_d d_mat(2, 2);
  d_mat << 10.0, 20.0, 30.0, 40.0;
  auto o3 = stan::math::make_partials_propagator(d_mat);

  EXPECT_EQ(2, sizeof(o3));

  matrix_v v_mat(2, 2);
  var v1 = var(0.0);
  var v2 = var(1.0);
  var v3 = var(2.0);
  var v4 = var(3.0);
  v_mat << v1, v2, v3, v4;

  std::vector<var> v_stdvec;
  v_stdvec.push_back(v1);
  v_stdvec.push_back(v2);
  v_stdvec.push_back(v3);
  v_stdvec.push_back(v4);

  auto o4 = stan::math::make_partials_propagator(v_mat);
  stan::math::edge<0>(o4).partials_ += d_mat;
  stan::math::edge<0>(o4).partials_vec_[1] += d_mat;
  // Should affect the same vars as the call above
  stan::math::edge<0>(o4).partials_vec_[27] += d_mat;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(30.0, grad[0]);
  EXPECT_FLOAT_EQ(60.0, grad[1]);
  EXPECT_FLOAT_EQ(90.0, grad[2]);
  EXPECT_FLOAT_EQ(120.0, grad[3]);
}

TEST(AgradPartialsVari, PartialsPropagatorMatMultivar) {
  using stan::math::make_partials_propagator;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::var;

  matrix_d d_mat(2, 2);
  d_mat << 10.0, 20.0, 30.0, 40.0;
  std::vector<matrix_d> d_mat_vec;
  d_mat_vec.push_back(d_mat);
  auto o3 = stan::math::make_partials_propagator(d_mat_vec);

  EXPECT_EQ(2, sizeof(o3));

  matrix_v v_mat1(2, 2);
  var v1 = var(0.0);
  var v2 = var(1.0);
  var v3 = var(2.0);
  var v4 = var(3.0);
  v_mat1 << v1, v2, v3, v4;

  matrix_v v_mat2(2, 2);
  var v5 = var(0.1);
  var v6 = var(1.1);
  var v7 = var(2.1);
  var v8 = var(3.1);
  v_mat2 << v5, v6, v7, v8;

  std::vector<matrix_v> v_mat_vec;
  v_mat_vec.push_back(v_mat1);
  v_mat_vec.push_back(v_mat2);

  std::vector<var> v_stdvec;
  v_stdvec.push_back(v1);
  v_stdvec.push_back(v2);
  v_stdvec.push_back(v3);
  v_stdvec.push_back(v4);
  v_stdvec.push_back(v5);
  v_stdvec.push_back(v6);
  v_stdvec.push_back(v7);
  v_stdvec.push_back(v8);

  auto o4 = stan::math::make_partials_propagator(v_mat_vec);
  stan::math::edge<0>(o4).partials_vec_[0] += d_mat;
  // Should NOT affect the same vars as the call above
  stan::math::edge<0>(o4).partials_vec_[1] += d_mat;
  stan::math::edge<0>(o4).partials_vec_[1] += d_mat;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
  EXPECT_FLOAT_EQ(20.0, grad[1]);
  EXPECT_FLOAT_EQ(30.0, grad[2]);
  EXPECT_FLOAT_EQ(40.0, grad[3]);
  EXPECT_FLOAT_EQ(20.0, grad[4]);
  EXPECT_FLOAT_EQ(40.0, grad[5]);
  EXPECT_FLOAT_EQ(60.0, grad[6]);
  EXPECT_FLOAT_EQ(80.0, grad[7]);
}

TEST(AgradPartialsVari, PartialsPropagatorMultivar) {
  using stan::math::make_partials_propagator;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  std::vector<vector_d> d_vec_vec(2);
  vector_d d_vec1(2);
  d_vec1 << 10.0, 20.0;
  vector_d d_vec2(2);
  d_vec2 << 30.0, 40.0;
  d_vec_vec.push_back(d_vec1);
  d_vec_vec.push_back(d_vec2);
  auto o3 = stan::math::make_partials_propagator(d_vec_vec);

  EXPECT_EQ(2, sizeof(o3));

  vector_v v_vec1(2);
  var v1 = var(0.0);
  var v2 = var(1.0);
  v_vec1 << v1, v2;
  vector_v v_vec2(2);
  var v3 = var(2.0);
  var v4 = var(3.0);
  v_vec2 << v3, v4;
  std::vector<vector_v> v_vec;
  v_vec.push_back(v_vec1);
  v_vec.push_back(v_vec2);

  std::vector<var> v_stdvec;
  v_stdvec.push_back(v1);
  v_stdvec.push_back(v2);
  v_stdvec.push_back(v3);
  v_stdvec.push_back(v4);

  auto o4 = stan::math::make_partials_propagator(v_vec);
  stan::math::edge<0>(o4).partials_vec_[0] += d_vec1;
  stan::math::edge<0>(o4).partials_vec_[1] += d_vec2;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
  EXPECT_FLOAT_EQ(20.0, grad[1]);
  EXPECT_FLOAT_EQ(30.0, grad[2]);
  EXPECT_FLOAT_EQ(40.0, grad[3]);
}

// XXX Test mixed - partials_propagator<std::vector<matrix_v>,
//                                        vector_d, vector_v>
TEST(AgradPartialsVari, PartialsPropagatorMultivarMixed) {
  using stan::math::make_partials_propagator;
  using stan::math::matrix_v;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  std::vector<vector_d> d_vec_vec(2);
  vector_d d_vec1(2);
  d_vec1 << 10.0, 20.0;
  vector_d d_vec2(2);
  d_vec2 << 30.0, 40.0;
  d_vec_vec.push_back(d_vec1);
  d_vec_vec.push_back(d_vec2);

  vector_v v_vec1(2);
  var v1 = var(0.0);
  var v2 = var(1.0);
  v_vec1 << v1, v2;
  vector_v v_vec2(2);
  var v3 = var(2.0);
  var v4 = var(3.0);
  v_vec2 << v3, v4;
  std::vector<vector_v> v_vec;
  v_vec.push_back(v_vec1);

  std::vector<var> v_stdvec;
  v_stdvec.push_back(v1);
  v_stdvec.push_back(v2);
  v_stdvec.push_back(v3);
  v_stdvec.push_back(v4);

  auto o4 = make_partials_propagator(v_vec, d_vec_vec, v_vec2);
  stan::math::edge<0>(o4).partials_vec_[0] += d_vec1;
  stan::math::edge<2>(o4).partials_vec_[0] += d_vec2;

  // 2 partials stdvecs, 4 pointers to edges, 2 pointers to operands
  // vecs
  EXPECT_EQ(112, sizeof(o4));

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
  EXPECT_FLOAT_EQ(20.0, grad[1]);
  EXPECT_FLOAT_EQ(30.0, grad[2]);
  EXPECT_FLOAT_EQ(40.0, grad[3]);

  // when given vector_d in place of vector_v all expressions must
  // still compile
  auto o5 = make_partials_propagator(v_vec, d_vec_vec, d_vec2);
  stan::math::edge<0>(o5).partials_vec_[0] += d_vec1;
  if (false) {
    // the test here is to make things compile as this pattern to
    // if-out things when terms are const is used in our functions
    stan::math::edge<2>(o5).partials_vec_[0] += vector_d();
    stan::math::edge<2>(o5).partials_vec_[0] -= vector_d();
    stan::math::edge<2>(o5).partials_vec_[0](0) = 0;
  }

  // the same needs to work for the nested case
  auto o6 = make_partials_propagator(d_vec_vec, d_vec_vec, v_vec2);
  if (false) {
    // the test here is to make things compile as this pattern to
    // if-out things when terms are const is used in our functions
    stan::math::edge<0>(o6).partials_vec_[0] += d_vec1;
  }
  stan::math::edge<2>(o6).partials_vec_[0] += d_vec2;
}

TEST(AgradPartialsVari, PartialsPropagatorVarValueMat) {
  using stan::math::make_partials_propagator;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::var;

  Eigen::MatrixXd a(2, 2);
  a << 10.0, 20.0, 30.0, 40.0;
  stan::math::var_value<Eigen::MatrixXd> av(a);

  auto ops = stan::math::make_partials_propagator(av);

  stan::math::edge<0>(ops).partials_ = Eigen::MatrixXd::Constant(2, 2, -2);
  var lp = ops.build(1);
  (2 * lp).grad();
  EXPECT_MATRIX_EQ(av.adj(), Eigen::MatrixXd::Constant(2, 2, -4))
}

TEST(AgradPartialsVari, PartialsPropagatorStdVectorVarValueMat) {
  using stan::math::make_partials_propagator;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::var;

  Eigen::MatrixXd a(2, 2);
  a << 10.0, 20.0, 30.0, 40.0;
  std::vector<stan::math::var_value<Eigen::MatrixXd>> av{a, a};

  auto ops = make_partials_propagator(av);

  stan::math::edge<0>(ops).partials_vec_[0]
      = Eigen::MatrixXd::Constant(2, 2, -2);
  stan::math::edge<0>(ops).partials_vec_[1]
      = Eigen::MatrixXd::Constant(2, 2, -3);
  var lp = ops.build(1);
  (2 * lp).grad();
  EXPECT_MATRIX_EQ(av[0].adj(), Eigen::MatrixXd::Constant(2, 2, -4));
  EXPECT_MATRIX_EQ(av[1].adj(), Eigen::MatrixXd::Constant(2, 2, -6));
}
