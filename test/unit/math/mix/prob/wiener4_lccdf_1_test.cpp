#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_y_a_t0, wiener4_lccdf) {
auto f_y_a_t0 = [](const auto& w, const auto& v) {
     return [&w, &v](const auto& y, const auto& a, const auto& t0) {
          return stan::math::wiener4_lccdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_y_a_t0(w, v), y, a, t0);
    }



TEST(mathMixScalFun_y_a_w, wiener4_lccdf) {
auto f_y_a_w = [](const auto& t0, const auto& v) {
     return [&t0, &v](const auto& y, const auto& a, const auto& w) {
          return stan::math::wiener4_lccdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_y_a_w(t0, v), y, a, w);
    }



TEST(mathMixScalFun_y_a_v, wiener4_lccdf) {
auto f_y_a_v = [](const auto& t0, const auto& w) {
     return [&t0, &w](const auto& y, const auto& a, const auto& v) {
          return stan::math::wiener4_lccdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_y_a_v(t0, w), y, a, v);
    }



TEST(mathMixScalFun_y_t0_w, wiener4_lccdf) {
auto f_y_t0_w = [](const auto& a, const auto& v) {
     return [&a, &v](const auto& y, const auto& t0, const auto& w) {
          return stan::math::wiener4_lccdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_y_t0_w(a, v), y, t0, w);
    }



TEST(mathMixScalFun_y_t0_v, wiener4_lccdf) {
auto f_y_t0_v = [](const auto& a, const auto& w) {
     return [&a, &w](const auto& y, const auto& t0, const auto& v) {
          return stan::math::wiener4_lccdf(y, a, t0, w, v, 1e-4);
    };
  };
 
      double y = 1.0;
      double a = 2.5;
      double t0 = 0.2;
      double w = 0.5;
      double v = 2.0;
      
      stan::test::expect_ad(f_y_t0_v(a, w), y, t0, v);
    }
