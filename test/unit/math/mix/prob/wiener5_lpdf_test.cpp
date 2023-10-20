#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>


TEST(mathMixDouble, wiener5_lpdf) { 
	double y = 1.0;
	double a = 2.0;
	double t0 = 0.2;
	double w = 0.5;
	double v = 1.5;
	double sv = 0.2;
	stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
}


TEST(mathMixVar, wiener5_lpdf) { 
	using stan::math::var;
	var y = 1.0;
	var a = 2.0;
	var t0 = 0.2;
	var w = 0.5;
	var v = 1.5;
	var sv = 0.2;
	stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
}


TEST(mathMixFVar, wiener5_lpdf) {
	using stan::math::fvar;
	using stan::math::var;
	fvar<var> y = 1.0;
	fvar<var> a = 2.0;
	fvar<var> t0 = 0.2;
	fvar<var> w = 0.5;
	fvar<var> v = 1.5;
	fvar<var> sv = 0.2;	
	double error = 1e-4;
	stan::math::wiener5_lpdf(y, a, t0, w, v, sv, error);
}



TEST(mathMixDouble1Fun, wiener5_lpdf) { 
	using stan::math::var;
	double y = 1.0;
	double a = 2.0;
	double t0 = 0.2;
	double w = 0.5;
	double v = 1.5;
	double sv = 0.2;
	auto f1 = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& v) {
	  return [&y, &a, &t0, &w, &v](const auto& sv) {
		return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
	  };
	};
	stan::test::expect_ad(f1(y, a, t0, w, v), sv);
}
