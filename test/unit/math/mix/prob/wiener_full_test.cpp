#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixDouble, wiener_full_lpdf) { 
	double y = 1.0;
	double a = 2.0;
	double t0 = 0.2;
	double w = 0.5;
	double v = 1.5;
	double sv = 0.2;
	double sw = 0.2; 
	double st0 = 0.2;
	stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0, 1e-4);
}


TEST(mathMixVar, wiener_full_lpdf) { 
	using stan::math::var;
	var y = 1.0;
	var a = 2.0;
	var t0 = 0.2;
	var w = 0.5;
	var v = 1.5;
	var sv = 0.2;
	var sw = 0; 
	var st0 = 0.2;
	stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
}


TEST(mathMixFVar, wiener_full_lpdf) { 
	using stan::math::var;
	using stan::math::fvar;
	fvar<var> y = 1.0;
	fvar<var> a = 2.0;
	fvar<var> t0 = 0.2;
	fvar<var> w = 0.5;
	fvar<var> v = 1.5;
	fvar<var> sv = 0.2;
	fvar<var> sw = 0; 
	fvar<var> st0 = 0.2;
	stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
}


TEST(mathMixScalFun, wiener_full_lpdf) { // ... runs within a few minutes when 1 parameter is tested ... TODO: have to add tests for the other 5 parameters
	
	auto f_sw = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& v, const auto& sv, const auto& st0) {
	  return [&y, &a, &t0, &w, &v, &sv, &st0](const auto& sw) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
	auto f_st0 = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& v, const auto& sv, const auto& sw) {
	  return [&y, &a, &t0, &w, &v, &sv, &sw](const auto& st0) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
	double y = 0.1;
	double a = 2.0;
	double t0 = 0.2;
	double w = 0.5;
	double v = 2.0;
	double sv = 0.2;
	double sw = 0.1;
	double st0 = 0.3;

	stan::test::expect_ad(f_st0(y, a, t0, w, v, sv, st0), sw);
	stan::test::expect_ad(f_st0(y, a, t0, w, v, sv, sw), st0);

}

/*
TEST(mathMixVecFun, wiener_full_lpdf) {
	auto f1 = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& v) {
	  return [&y, &a, &t0, &w, &v](const auto& sv, const auto& sw, const auto& st0) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
/*	auto f2 = [](const auto& y, const auto& a, const auto& sv, const auto& sw, const auto& st0) {
	  return [&y, &a, &sv, &sw, &st0](const auto& t0, const auto& w, const auto& v) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
	auto f3 = [](const auto& w, const auto& v, const auto& sv, const auto& sw, const auto& st0) {
	  return [&w, &v, &sv, &sw, &st0](const auto& y, const auto& a, const auto& t0) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};

	Eigen::VectorXd y(2);
	y << 1.0, 1.5;
	Eigen::VectorXd a(2);
	a << 2.0, 2.5;
	Eigen::VectorXd t0(2);
	t0 << 0.2, 0.1;
	Eigen::VectorXd w(2);
	w << 0.5, 0.45;
	Eigen::VectorXd v(2);
	v << 2.0, 1.0;
	Eigen::VectorXd sv(2);
	sv << 0.2, 0.1;
	Eigen::VectorXd sw(2);
	sw << 0.1, 0.2;
	Eigen::VectorXd st0(2);
	st0 << 0.3, 0.2;

	stan::test::expect_ad(f1(y, a, t0, w, v), sv, sw, st0);
//	stan::test::expect_ad(f2(y, a, sv, sw, st0), t0, w, v);
//	stan::test::expect_ad(f3(w, v, sv, sw, st0), y, a, t0);


}
*/
