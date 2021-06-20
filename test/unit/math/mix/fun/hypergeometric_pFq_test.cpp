#include <test/unit/math/test_ad.hpp>
#include <limits>


TEST(mathMixScalFun, hyper_pfq) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::hypergeometric_pFq(x1, x2, x3);
  };

  Eigen::VectorXd p(2);
  p << 4.0, 2.0;

  Eigen::VectorXd q(2);
  q << 6.0, 3.0;

  double z = 4.0;

  stan::test::expect_ad(f, p, q, z);
}


TEST(primScalFun, grad_2F2_ffv) {
  using stan::math::vector_ffv;
  using stan::math::vector_d;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::hypergeometric_pFq;

  vector_ffv fd_p(2);
  vector_d d_p(2);
  fd_p.val() << 4, 2;
  d_p << 4, 2;
  fd_p.val().d() << 1, 1;

  vector_ffv fd_q(2);
  vector_d d_q(2);
  fd_q.val() << 6, 3;
  d_q << 6, 3;
  fd_q.val().d() << 1, 1;

  fvar<fvar<var>> ffd_z;
  ffd_z.val_.val_= 4;
  ffd_z.val_.d_ = 1;
  double d_z = 4;

  // fvar, fvar, fvar
  fvar<fvar<var>> result = hypergeometric_pFq(d_p, fd_q, d_z);

  double adj = 3.924636646666071 + 6.897245961898751
                 -2.775051002566842 - 4.980095849781222
                 + 4.916522138006060;
  
  //EXPECT_FLOAT_EQ(result.val_.d_, adj);
}
