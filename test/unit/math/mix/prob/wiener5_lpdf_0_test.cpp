#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixFVarFvarVar, wiener5_lpdf) {
  using stan::math::fvar;
  using stan::math::var;
  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  fvar<fvar<var>> w(fvar<var>(0.5, 1.0), fvar<var>(1.0, 0));
  double v = 1.5;
  double sv = 0.2;
  double error = 1e-4;
  fvar<fvar<var>> res = stan::math::wiener5_lpdf(y, a, t0, w, v, sv, error);
  grad(res.d_.d_.vi_);
  std::cout << "res" << std::endl;
  std::cout << res.val().val().val() << std::endl;
  std::cout << res.val().d().val() << std::endl;
  std::cout << res.d().val().val() << std::endl;
  std::cout << res.d().d().val() << std::endl;

#ifdef NOPENO
std::cout << "\ny" << std::endl;
  std::cout << y.val().val().adj() << std::endl;
  std::cout << y.val().d().adj() << std::endl;
  std::cout << y.d().val().adj() << std::endl;
  std::cout << y.d().d().adj() << std::endl;
std::cout << "\na" << std::endl;
  std::cout << a.val().val().adj() << std::endl;
  std::cout << a.val().d().adj() << std::endl;
  std::cout << a.d().val().adj() << std::endl;
  std::cout << a.d().d().adj() << std::endl;

std::cout << "\nt0" << std::endl;
  std::cout << t0.val().val().adj() << std::endl;
  std::cout << t0.val().d().adj() << std::endl;
  std::cout << t0.d().val().adj() << std::endl;
  std::cout << t0.d().d().adj() << std::endl;
#endif

std::cout << "\nw" << std::endl;
  std::cout << w.val().val().adj() << std::endl;
  std::cout << w.val().d().adj() << std::endl;
  std::cout << w.d().val().adj() << std::endl;
  std::cout << w.d().d().adj() << std::endl;
#ifdef NOPENO
std::cout << "\nv" << std::endl;
  std::cout << v.val().val().adj() << std::endl;
  std::cout << v.val().d().adj() << std::endl;
  std::cout << v.d().val().adj() << std::endl;
  std::cout << v.d().d().adj() << std::endl;

std::cout << "\nsv" << std::endl;
  std::cout << sv.val().val().adj() << std::endl;
  std::cout << sv.val().d().adj() << std::endl;
  std::cout << sv.d().val().adj() << std::endl;
  std::cout << sv.d().d().adj() << std::endl;
#endif
}

