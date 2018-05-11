#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradPartialsVari, OperandsAndPartialsUniMixMat) {
  using stan::math::fvar;
  using stan::math::operands_and_partials;
  using stan::math::var;

  std::vector<var> val_dxs;
  val_dxs.push_back(1.0);
  val_dxs.push_back(2.0);
  val_dxs.push_back(3.0);
  val_dxs.push_back(4.0);

  Eigen::Matrix<fvar<var>, -1, -1> m1(2, 2);
  // Set d_ to 1 for one variable we care about;
  m1 << fvar<var>(val_dxs[0], 1.0), fvar<var>(val_dxs[1], 0.0),
      fvar<var>(val_dxs[2], 0.0), fvar<var>(val_dxs[3], 0.0);

  Eigen::Matrix<var, -1, -1> dxm1(2, 2);
  dxm1 << 4.0, 5.0, 6.0, 7.0;
  std::vector<var> d_dxs;
  d_dxs.push_back(dxm1(0));
  d_dxs.push_back(dxm1(1));
  d_dxs.push_back(dxm1(2));
  d_dxs.push_back(dxm1(3));

  operands_and_partials<Eigen::Matrix<fvar<var>, -1, -1>> ops_partials(m1);

  // Do normal math on the fvar<var>, do addition to partials,
  ops_partials.edge1_.partials_vec_[0] += dxm1;

  std::vector<double> gradient;
  fvar<var> v = ops_partials.build(10.0);

  // call grad(f.d_) afterwards
  v.d_.grad(d_dxs, gradient);

  EXPECT_FLOAT_EQ(10.0, v.val().val());
  EXPECT_FLOAT_EQ(4.0, v.d_.val());
  EXPECT_FLOAT_EQ(1.0, gradient[0]);
  EXPECT_FLOAT_EQ(0.0, gradient[1]);
  EXPECT_FLOAT_EQ(0.0, gradient[2]);
  EXPECT_FLOAT_EQ(0.0, gradient[3]);
}

TEST(AgradPartialsVari, OperandsAndPartialsUniMixMat_dbl) {
  using stan::is_constant_struct;
  using stan::math::fvar;
  using stan::math::operands_and_partials;
  using stan::math::var;

  std::vector<var> val_dxs;
  val_dxs.push_back(1.0);
  val_dxs.push_back(2.0);
  val_dxs.push_back(3.0);
  val_dxs.push_back(4.0);

  Eigen::Matrix<fvar<var>, -1, -1> m1(2, 2);
  m1 << fvar<var>(val_dxs[0], 1.0), fvar<var>(val_dxs[1], 0.0),
      fvar<var>(val_dxs[2], 0.0), fvar<var>(val_dxs[3], 0.0);

  Eigen::Matrix<double, -1, -1> m2(2, 2);
  m2 << 2.0, 2.0, 2.0, 2.0;

  Eigen::Matrix<var, -1, -1> dxm1(2, 2);
  dxm1 << 4.0, 5.0, 6.0, 7.0;
  std::vector<var> d_dxs;
  d_dxs.push_back(dxm1(0));
  d_dxs.push_back(dxm1(1));
  d_dxs.push_back(dxm1(2));
  d_dxs.push_back(dxm1(3));

  operands_and_partials<Eigen::Matrix<fvar<var>, -1, -1>,
                        Eigen::Matrix<double, -1, -1>>
      ops_partials(m1, m2);

  // Do normal math on the fvar<var>, do addition to partials,
  ops_partials.edge1_.partials_vec_[0] += dxm1;

  std::vector<double> gradient;
  fvar<var> v = ops_partials.build(10.0);

  // call grad(f.d_) afterwards
  v.d_.grad(d_dxs, gradient);

  EXPECT_FLOAT_EQ(10.0, v.val().val());
  EXPECT_FLOAT_EQ(4.0, v.d_.val());
  EXPECT_FLOAT_EQ(1.0, gradient[0]);
  EXPECT_FLOAT_EQ(0.0, gradient[1]);
  EXPECT_FLOAT_EQ(0.0, gradient[2]);
  EXPECT_FLOAT_EQ(0.0, gradient[3]);
}

TEST(AgradPartialsVari, OperandsAndPartialsMultiMix) {
  using stan::math::fvar;
  using stan::math::operands_and_partials;
  using stan::math::var;

  typedef Eigen::Matrix<fvar<var>, -1, -1> uni_mat_t;
  std::vector<var> val_dxs;
  val_dxs.push_back(1.0);
  val_dxs.push_back(2.0);
  val_dxs.push_back(3.0);
  val_dxs.push_back(4.0);
  val_dxs.push_back(5.0);
  val_dxs.push_back(6.0);
  val_dxs.push_back(7.0);
  val_dxs.push_back(8.0);

  uni_mat_t m1(2, 2);
  // Set d_ to 1 for one variable we care about;
  m1(0, 0) = val_dxs[0];
  m1(0, 0).d_ = 1;
  m1(0, 1) = val_dxs[1];
  m1(1, 0) = val_dxs[2];
  m1(1, 1) = val_dxs[3];

  uni_mat_t m2(2, 2);
  m2(0, 0) = val_dxs[4];
  m2(0, 1) = val_dxs[5];
  m2(1, 0) = val_dxs[6];
  m2(1, 1) = val_dxs[7];

  std::vector<uni_mat_t> multi_mat;
  multi_mat.push_back(m1);
  multi_mat.push_back(m2);

  Eigen::Matrix<var, -1, -1> dxm1(2, 2);
  dxm1 << 4.0, 5.0, 6.0, 7.0;
  std::vector<var> d_dxs;
  d_dxs.push_back(dxm1(0));
  d_dxs.push_back(dxm1(1));
  d_dxs.push_back(dxm1(2));
  d_dxs.push_back(dxm1(3));

  operands_and_partials<std::vector<uni_mat_t>> ops_partials(multi_mat);

  // Do normal math on the fvar<var>, do addition to partials,
  ops_partials.edge1_.partials_vec_[0] += dxm1;
  ops_partials.edge1_.partials_vec_[1] += 2 * dxm1;

  std::vector<double> gradient;
  fvar<var> v = ops_partials.build(10.0);

  // call grad(f.d_) afterwards
  v.d_.grad(d_dxs, gradient);

  EXPECT_FLOAT_EQ(10.0, v.val().val());
  EXPECT_FLOAT_EQ(4.0, v.d_.val());
  EXPECT_FLOAT_EQ(1.0, gradient[0]);
  EXPECT_FLOAT_EQ(0.0, gradient[1]);
  EXPECT_FLOAT_EQ(0.0, gradient[2]);
  EXPECT_FLOAT_EQ(0.0, gradient[3]);
}

TEST(AgradPartialsVari, OperandsAndPartialsMultiMix_dbl) {
  using stan::math::fvar;
  using stan::math::operands_and_partials;
  using stan::math::var;

  typedef Eigen::Matrix<fvar<var>, -1, -1> uni_mat_t;
  std::vector<var> val_dxs;
  val_dxs.push_back(1.0);
  val_dxs.push_back(2.0);
  val_dxs.push_back(3.0);
  val_dxs.push_back(4.0);
  val_dxs.push_back(5.0);
  val_dxs.push_back(6.0);
  val_dxs.push_back(7.0);
  val_dxs.push_back(8.0);

  uni_mat_t m1(2, 2);
  // Set d_ to 1 for one variable we care about;
  m1(0, 0) = val_dxs[0];
  m1(0, 0).d_ = 1;
  m1(0, 1) = val_dxs[1];
  m1(1, 0) = val_dxs[2];
  m1(1, 1) = val_dxs[3];

  uni_mat_t m2(2, 2);
  m2(0, 0) = val_dxs[4];
  m2(0, 1) = val_dxs[5];
  m2(1, 0) = val_dxs[6];
  m2(1, 1) = val_dxs[7];

  std::vector<uni_mat_t> multi_mat;
  multi_mat.push_back(m1);
  multi_mat.push_back(m2);

  Eigen::Matrix<double, -1, -1> m3(2, 2);
  m3 << 2.0, 2.0, 2.0, 2.0;

  Eigen::Matrix<var, -1, -1> dxm1(2, 2);
  dxm1 << 4.0, 5.0, 6.0, 7.0;
  std::vector<var> d_dxs;
  d_dxs.push_back(dxm1(0));
  d_dxs.push_back(dxm1(1));
  d_dxs.push_back(dxm1(2));
  d_dxs.push_back(dxm1(3));

  operands_and_partials<std::vector<uni_mat_t>, Eigen::Matrix<double, -1, -1>>
      ops_partials(multi_mat, m3);

  // Do normal math on the fvar<var>, do addition to partials,
  ops_partials.edge1_.partials_vec_[0] += dxm1;
  ops_partials.edge1_.partials_vec_[1] += 2 * dxm1;

  std::vector<double> gradient;
  fvar<var> v = ops_partials.build(10.0);

  // call grad(f.d_) afterwards
  v.d_.grad(d_dxs, gradient);

  EXPECT_FLOAT_EQ(10.0, v.val().val());
  EXPECT_FLOAT_EQ(4.0, v.d_.val());
  EXPECT_FLOAT_EQ(1.0, gradient[0]);
  EXPECT_FLOAT_EQ(0.0, gradient[1]);
  EXPECT_FLOAT_EQ(0.0, gradient[2]);
  EXPECT_FLOAT_EQ(0.0, gradient[3]);
}

TEST(AgradPartialsVari, OperandsAndPartialsMultiStdMix) {
  using stan::math::fvar;
  using stan::math::operands_and_partials;
  using stan::math::var;

  typedef std::vector<fvar<var>> uni_std_t;
  std::vector<var> val_dxs;
  val_dxs.push_back(1.0);
  val_dxs.push_back(2.0);
  val_dxs.push_back(3.0);
  val_dxs.push_back(4.0);
  val_dxs.push_back(5.0);
  val_dxs.push_back(6.0);
  val_dxs.push_back(7.0);
  val_dxs.push_back(8.0);

  uni_std_t m1(4);
  // Set d_ to 1 for one variable we care about;
  m1[0] = val_dxs[0];
  m1[0].d_ = 1;
  m1[1] = val_dxs[1];
  m1[2] = val_dxs[2];
  m1[3] = val_dxs[3];

  uni_std_t m2(4);
  m2[0] = val_dxs[4];
  m2[1] = val_dxs[5];
  m2[2] = val_dxs[6];
  m2[3] = val_dxs[7];

  std::vector<uni_std_t> multi_std;
  multi_std.push_back(m1);
  multi_std.push_back(m2);

  std::vector<var> dxm1(4);
  dxm1[0] = 4.0;
  dxm1[1] = 5.0;
  dxm1[2] = 6.0;
  dxm1[3] = 7.0;

  std::vector<var> d_dxs;
  d_dxs.push_back(dxm1[0]);
  d_dxs.push_back(dxm1[1]);
  d_dxs.push_back(dxm1[2]);
  d_dxs.push_back(dxm1[3]);

  operands_and_partials<std::vector<uni_std_t>> ops_partials(multi_std);

  // Do normal math on the fvar<var>, do addition to partials,
  ops_partials.edge1_.partials_vec_[0] = dxm1;
  ops_partials.edge1_.partials_vec_[1] = dxm1;

  std::vector<double> gradient;
  fvar<var> v = ops_partials.build(10.0);

  // call grad(f.d_) afterwards
  v.d_.grad(d_dxs, gradient);

  EXPECT_FLOAT_EQ(10.0, v.val().val());
  EXPECT_FLOAT_EQ(4.0, v.d_.val());
  EXPECT_FLOAT_EQ(1.0, gradient[0]);
  EXPECT_FLOAT_EQ(0.0, gradient[1]);
  EXPECT_FLOAT_EQ(0.0, gradient[2]);
  EXPECT_FLOAT_EQ(0.0, gradient[3]);
}

TEST(AgradPartialsVari, OperandsAndPartialsMultiStdMix_dbl) {
  using stan::math::fvar;
  using stan::math::operands_and_partials;
  using stan::math::var;

  typedef std::vector<fvar<var>> uni_std_t;
  std::vector<var> val_dxs;
  val_dxs.push_back(1.0);
  val_dxs.push_back(2.0);
  val_dxs.push_back(3.0);
  val_dxs.push_back(4.0);
  val_dxs.push_back(5.0);
  val_dxs.push_back(6.0);
  val_dxs.push_back(7.0);
  val_dxs.push_back(8.0);

  uni_std_t m1(4);
  // Set d_ to 1 for one variable we care about;
  m1[0] = val_dxs[0];
  m1[0].d_ = 1;
  m1[1] = val_dxs[1];
  m1[2] = val_dxs[2];
  m1[3] = val_dxs[3];

  uni_std_t m2(4);
  m2[0] = val_dxs[4];
  m2[1] = val_dxs[5];
  m2[2] = val_dxs[6];
  m2[3] = val_dxs[7];

  std::vector<uni_std_t> multi_std;
  multi_std.push_back(m1);
  multi_std.push_back(m2);

  Eigen::Matrix<double, -1, -1> m3(2, 2);
  m3 << 2.0, 2.0, 2.0, 2.0;

  std::vector<var> dxm1(4);
  dxm1[0] = 4.0;
  dxm1[1] = 5.0;
  dxm1[2] = 6.0;
  dxm1[3] = 7.0;

  std::vector<var> d_dxs;
  d_dxs.push_back(dxm1[0]);
  d_dxs.push_back(dxm1[1]);
  d_dxs.push_back(dxm1[2]);
  d_dxs.push_back(dxm1[3]);

  operands_and_partials<std::vector<uni_std_t>, Eigen::Matrix<double, -1, -1>>
      ops_partials(multi_std, m3);

  // Do normal math on the fvar<var>, do addition to partials,
  ops_partials.edge1_.partials_vec_[0] = dxm1;
  ops_partials.edge1_.partials_vec_[1] = dxm1;

  std::vector<double> gradient;
  fvar<var> v = ops_partials.build(10.0);

  // call grad(f.d_) afterwards
  v.d_.grad(d_dxs, gradient);

  EXPECT_FLOAT_EQ(10.0, v.val().val());
  EXPECT_FLOAT_EQ(4.0, v.d_.val());
  EXPECT_FLOAT_EQ(1.0, gradient[0]);
  EXPECT_FLOAT_EQ(0.0, gradient[1]);
  EXPECT_FLOAT_EQ(0.0, gradient[2]);
  EXPECT_FLOAT_EQ(0.0, gradient[3]);
}

TEST(AgradPartialsVari, OperandsAndPartialsMultiMixInt) {
  using stan::math::operands_and_partials;

  typedef Eigen::Matrix<int, -1, -1> uni_mat_t;

  uni_mat_t m1(2, 2);
  // Set d_ to 1 for one variable we care about;
  m1(0, 0) = 1;
  m1(0, 1) = 0;
  m1(1, 0) = 0;
  m1(1, 1) = 1;

  uni_mat_t m2(2, 2);
  m2(0, 0) = 2;
  m2(0, 1) = 3;
  m2(1, 0) = 4;
  m2(1, 1) = 5;

  std::vector<uni_mat_t> multi_mat;
  multi_mat.push_back(m1);
  multi_mat.push_back(m2);

  operands_and_partials<std::vector<uni_mat_t>> ops_partials(multi_mat);

  double v = ops_partials.build(10.0);

  EXPECT_FLOAT_EQ(10.0, v);
}

TEST(AgradPartialsVari, OperandsAndPartialsMultiMixInt_dbl) {
  using stan::math::operands_and_partials;

  typedef Eigen::Matrix<int, -1, -1> uni_mat_t;

  uni_mat_t m1(2, 2);
  // Set d_ to 1 for one variable we care about;
  m1(0, 0) = 1;
  m1(0, 1) = 0;
  m1(1, 0) = 0;
  m1(1, 1) = 1;

  uni_mat_t m2(2, 2);
  m2(0, 0) = 2;
  m2(0, 1) = 3;
  m2(1, 0) = 4;
  m2(1, 1) = 5;

  std::vector<uni_mat_t> multi_mat;
  multi_mat.push_back(m1);
  multi_mat.push_back(m2);

  Eigen::Matrix<double, -1, -1> m3(2, 2);
  m3 << 2.0, 2.0, 2.0, 2.0;

  operands_and_partials<std::vector<uni_mat_t>, Eigen::Matrix<double, -1, -1>>
      ops_partials(multi_mat);

  double v = ops_partials.build(10.0);

  EXPECT_FLOAT_EQ(10.0, v);
}
