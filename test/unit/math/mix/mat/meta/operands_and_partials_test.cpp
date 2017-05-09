#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>

TEST(AgradPartialsVari, OperandsAndPartialsVec) {
  using stan::math::operands_and_partials;
  using stan::math::var;
  using stan::math::fvar;

  typedef fvar<var> scalar;
  typedef var dx;
  std::vector<dx> val_dxs;
  val_dxs.push_back(1.0);
  val_dxs.push_back(2.0);
  val_dxs.push_back(3.0);
  val_dxs.push_back(4.0);

  Eigen::Matrix<scalar, -1, -1> m1(2, 2);
  // Set d_ to 1 for one variable we care about;
  m1 << scalar(val_dxs[0], 1.0), scalar(val_dxs[1], 0.0),
    scalar(val_dxs[2], 0.0), scalar(val_dxs[3], 0.0);

  Eigen::Matrix<dx, -1, -1> dxm1(2, 2);
  dxm1 << 4.0, 5.0, 6.0, 7.0;
  std::vector<dx> d_dxs;
  d_dxs.push_back(dxm1(0));
  d_dxs.push_back(dxm1(1));
  d_dxs.push_back(dxm1(2));
  d_dxs.push_back(dxm1(3));

  operands_and_partials<Eigen::Matrix<scalar, -1, -1> > o1(m1);

  // Do normal math on the fvar<var>, do addition to partials,
  o1.edge1_.partials_vec_[0] += dxm1;

  std::vector<double> grad_grad;
  scalar v = o1.build(10.0);

  // call grad(f.d_) afterwards
  v.d_.grad(d_dxs, grad_grad);

  EXPECT_FLOAT_EQ(10.0, v.val().val());
  EXPECT_FLOAT_EQ(4.0, v.d_.val());
  EXPECT_FLOAT_EQ(1.0, grad_grad[0]);
  EXPECT_FLOAT_EQ(0.0, grad_grad[1]);
  EXPECT_FLOAT_EQ(0.0, grad_grad[2]);
  EXPECT_FLOAT_EQ(0.0, grad_grad[3]);
}
