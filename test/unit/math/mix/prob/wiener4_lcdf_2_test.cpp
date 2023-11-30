#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_y_w_v, wiener4_lcdf) {
auto f_y_w_v = [](const auto& a, const auto& t0) {
     return [&a, &t0](const auto& y, const auto& w, const auto& v) {
          return stan::math::wiener4_lcdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_y_w_v(a, t0), y, w, v);
    }



TEST(mathMixScalFun_a_t0_w, wiener4_lcdf) {
auto f_a_t0_w = [](const auto& y, const auto& v) {
     return [&y, &v](const auto& a, const auto& t0, const auto& w) {
          return stan::math::wiener4_lcdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_a_t0_w(y, v), a, t0, w);
    }



TEST(mathMixScalFun_a_t0_v, wiener4_lcdf) {
auto f_a_t0_v = [](const auto& y, const auto& w) {
     return [&y, &w](const auto& a, const auto& t0, const auto& v) {
          return stan::math::wiener4_lcdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_a_t0_v(y, w), a, t0, v);
    }



TEST(mathMixScalFun_a_w_v, wiener4_lcdf) {
auto f_a_w_v = [](const auto& y, const auto& t0) {
     return [&y, &t0](const auto& a, const auto& w, const auto& v) {
          return stan::math::wiener4_lcdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_a_w_v(y, t0), a, w, v);
    }



TEST(mathMixScalFun_t0_w_v, wiener4_lcdf) {
auto f_t0_w_v = [](const auto& y, const auto& a) {
     return [&y, &a](const auto& t0, const auto& w, const auto& v) {
          return stan::math::wiener4_lcdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_t0_w_v(y, a), t0, w, v);
    }
