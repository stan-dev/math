#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradPartialsVari, OperandsAndPartialsScal) {
  using stan::math::operands_and_partials;
  using stan::math::var;

  double d1;
  operands_and_partials<double> o3(d1);
  EXPECT_EQ(16, sizeof(o3));

  var v1 = var(0.0);

  std::vector<var> v_stdvec;
  v_stdvec.push_back(v1);

  operands_and_partials<var> o4(v1);
  o4.edge1_.partials_[0] += 10.0;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
}

TEST(AgradPartialsVari, OperandsAndPartialsVec) {
  using stan::math::operands_and_partials;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d_vec(4);
  operands_and_partials<vector_d> o3(d_vec);
  EXPECT_EQ(16, sizeof(o3));

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

  operands_and_partials<vector_v> o4(v_vec);
  o4.edge1_.partials_[0] += 10.0;
  o4.edge1_.partials_[1] += 20.0;
  o4.edge1_.partials_[2] += 30.0;
  o4.edge1_.partials_[3] += 40.0;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
  EXPECT_FLOAT_EQ(20.0, grad[1]);
  EXPECT_FLOAT_EQ(30.0, grad[2]);
  EXPECT_FLOAT_EQ(40.0, grad[3]);
}

TEST(AgradPartialsVari, OperandsAndPartialsStdVec) {
  using stan::math::operands_and_partials;
  using stan::math::var;

  std::vector<double> d_vec(4);
  operands_and_partials<std::vector<double>> o3(d_vec);
  EXPECT_EQ(16, sizeof(o3));

  std::vector<var> v_vec;
  var v1 = var(0.0);
  var v2 = var(1.0);
  var v3 = var(2.0);
  var v4 = var(3.0);
  v_vec.push_back(v1);
  v_vec.push_back(v2);
  v_vec.push_back(v3);
  v_vec.push_back(v4);

  operands_and_partials<std::vector<var>> o4(v_vec);
  o4.edge1_.partials_[0] += 10.0;
  o4.edge1_.partials_[1] += 20.0;
  o4.edge1_.partials_[2] += 30.0;
  o4.edge1_.partials_[3] += 40.0;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_vec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
  EXPECT_FLOAT_EQ(20.0, grad[1]);
  EXPECT_FLOAT_EQ(30.0, grad[2]);
  EXPECT_FLOAT_EQ(40.0, grad[3]);
}

TEST(AgradPartialsVari, OperandsAndPartialsMat) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::operands_and_partials;
  using stan::math::var;

  matrix_d d_mat(2, 2);
  d_mat << 10.0, 20.0, 30.0, 40.0;
  operands_and_partials<matrix_d> o3(d_mat);

  EXPECT_EQ(16, sizeof(o3));

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

  operands_and_partials<matrix_v> o4(v_mat);
  o4.edge1_.partials_ += d_mat;
  o4.edge1_.partials_vec_[1] += d_mat;
  // Should affect the same vars as the call above
  o4.edge1_.partials_vec_[27] += d_mat;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(30.0, grad[0]);
  EXPECT_FLOAT_EQ(60.0, grad[1]);
  EXPECT_FLOAT_EQ(90.0, grad[2]);
  EXPECT_FLOAT_EQ(120.0, grad[3]);
}

TEST(AgradPartialsVari, OperandsAndPartialsMatMultivar) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::operands_and_partials;
  using stan::math::var;

  matrix_d d_mat(2, 2);
  d_mat << 10.0, 20.0, 30.0, 40.0;
  std::vector<matrix_d> d_mat_vec;
  d_mat_vec.push_back(d_mat);
  operands_and_partials<std::vector<matrix_d>> o3(d_mat_vec);

  EXPECT_EQ(16, sizeof(o3));

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

  operands_and_partials<std::vector<matrix_v>> o4(v_mat_vec);
  o4.edge1_.partials_vec_[0] += d_mat;
  // Should NOT affect the same vars as the call above
  o4.edge1_.partials_vec_[1] += d_mat;
  o4.edge1_.partials_vec_[1] += d_mat;

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

TEST(AgradPartialsVari, OperandsAndPartialsMultivar) {
  using stan::math::operands_and_partials;
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
  operands_and_partials<std::vector<vector_d>> o3(d_vec_vec);

  EXPECT_EQ(16, sizeof(o3));

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

  operands_and_partials<std::vector<vector_v>> o4(v_vec);
  o4.edge1_.partials_vec_[0] += d_vec1;
  o4.edge1_.partials_vec_[1] += d_vec2;

  std::vector<double> grad;
  var v = o4.build(10.0);
  v.grad(v_stdvec, grad);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
  EXPECT_FLOAT_EQ(20.0, grad[1]);
  EXPECT_FLOAT_EQ(30.0, grad[2]);
  EXPECT_FLOAT_EQ(40.0, grad[3]);
}

// XXX Test mixed - operands_and_partials<std::vector<matrix_v>,
//                                        vector_d, vector_v>
TEST(AgradPartialsVari, OperandsAndPartialsMultivarMixed) {
  using stan::math::matrix_v;
  using stan::math::operands_and_partials;
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

  operands_and_partials<std::vector<vector_v>, std::vector<vector_d>, vector_v>
      o4(v_vec, d_vec_vec, v_vec2);
  o4.edge1_.partials_vec_[0] += d_vec1;
  o4.edge3_.partials_vec_[0] += d_vec2;

  // 2 partials stdvecs, 4 pointers to edges, 2 pointers to operands
  // vecs
  EXPECT_EQ(128, sizeof(o4));

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
  operands_and_partials<std::vector<vector_v>, std::vector<vector_d>, vector_d>
      o5(v_vec, d_vec_vec, d_vec2);
  o5.edge1_.partials_vec_[0] += d_vec1;
  if (false) {
    // the test here is to make things compile as this pattern to
    // if-out things when terms are const is used in our functions
    o5.edge3_.partials_vec_[0] += vector_d();
    o5.edge3_.partials_vec_[0] -= vector_d();
    o5.edge3_.partials_vec_[0](0) = 0;
  }

  // the same needs to work for the nested case
  operands_and_partials<std::vector<vector_d>, std::vector<vector_d>, vector_v>
      o6(d_vec_vec, d_vec_vec, v_vec2);
  if (false) {
    // the test here is to make things compile as this pattern to
    // if-out things when terms are const is used in our functions
    o6.edge1_.partials_vec_[0] += d_vec1;
  }
  o6.edge3_.partials_vec_[0] += d_vec2;
}

TEST(AgradPartialsVari, OperandsAndPartialsVarValueMat) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::operands_and_partials;
  using stan::math::var;

  Eigen::MatrixXd a(2, 2);
  a << 10.0, 20.0, 30.0, 40.0;
  stan::math::var_value<Eigen::MatrixXd> av(a);

  operands_and_partials<stan::math::var_value<Eigen::MatrixXd>> ops(av);

  ops.edge1_.partials_ = Eigen::MatrixXd::Constant(2, 2, -2);
  var lp = ops.build(1);
  (2 * lp).grad();
  EXPECT_MATRIX_EQ(av.adj(), Eigen::MatrixXd::Constant(2, 2, -4))
}

TEST(AgradPartialsVari, OperandsAndPartialsStdVectorVarValueMat) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::operands_and_partials;
  using stan::math::var;

  Eigen::MatrixXd a(2, 2);
  a << 10.0, 20.0, 30.0, 40.0;
  std::vector<stan::math::var_value<Eigen::MatrixXd>> av{a, a};

  operands_and_partials<std::vector<stan::math::var_value<Eigen::MatrixXd>>>
      ops(av);

  ops.edge1_.partials_vec_[0] = Eigen::MatrixXd::Constant(2, 2, -2);
  ops.edge1_.partials_vec_[1] = Eigen::MatrixXd::Constant(2, 2, -3);
  var lp = ops.build(1);
  (2 * lp).grad();
  EXPECT_MATRIX_EQ(av[0].adj(), Eigen::MatrixXd::Constant(2, 2, -4));
  EXPECT_MATRIX_EQ(av[1].adj(), Eigen::MatrixXd::Constant(2, 2, -6));
}
