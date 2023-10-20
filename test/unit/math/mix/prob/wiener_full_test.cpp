#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

/*
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
*/

TEST(mathMixScalFun, wiener_full_lpdf) { 
	
/*	auto f_y = [](const auto& a, const auto& t0, const auto& w, const auto& v, const auto& sv, const auto& sw, const auto& st0) {
	  return [&a, &t0, &w, &v, &sv, &sw, &st0](const auto& y) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
	auto f_a = [](const auto& y, const auto& t0, const auto& w, const auto& v, const auto& sv, const auto& sw, const auto& st0) {
	  return [&y, &t0, &w, &v, &sv, &sw, &st0](const auto& a) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
	auto f_t0 = [](const auto& y, const auto& a, const auto& w, const auto& v, const auto& sv, const auto& sw, const auto& st0) {
	  return [&y, &a, &w, &v, &sv, &sw, &st0](const auto& t0) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
	auto f_w = [](const auto& y, const auto& a, const auto& t0, const auto& v, const auto& sv, const auto& sw, const auto& st0) {
	  return [&y, &a, &t0, &v, &sv, &sw, &st0](const auto& w) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};

	auto f_v = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& sv, const auto& sw, const auto& st0) {
	  return [&y, &a, &t0, &w, &sv, &sw, &st0](const auto& v) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
	auto f_sv = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& v, const auto& sw, const auto& st0) {
	  return [&y, &a, &t0, &w, &v, &sw, &st0](const auto& sv) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
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
	
	*/
	auto f1 = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& v) {
	  return [&y, &a, &t0, &w, &v](const auto& sv, const auto& sw, const auto& st0) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	/*
	auto f2 = [](const auto& y, const auto& a, const auto& sv, const auto& sw, const auto& st0) {
	  return [&y, &a, &sv, &sw, &st0](const auto& t0, const auto& w, const auto& v) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
	
	auto f3 = [](const auto& w, const auto& v, const auto& sv, const auto& sw, const auto& st0) {
	  return [&w, &v, &sv, &sw, &st0](const auto& y, const auto& a, const auto& t0) {
		return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
	  };
	};
*/

	double y = 0.1;
	double a = 2.0;
	double t0 = 0.2;
	double w = 0.5;
	double v = 2.0;
	double sv = 0.2;
	double sw = 0.1;
	double st0 = 0.3;
	

	
	stan::test::expect_ad(f1(y, a, t0, w, v), sv, sw, st0);
//	stan::test::expect_ad(f2(y, a, sv, sw, st0), t0, w, v);
//	stan::test::expect_ad(f3(w, v, sv, sw, st0), y, a, t0);
	

	
	
//	stan::test::expect_ad(f_y(a, t0, w, v, sv, sw, st0), y);
//	stan::test::expect_ad(f_a(y, t0, w, v, sv, sw, st0), a);
//	stan::test::expect_ad(f_t0(y, a, w, v, sv, sw, st0), t0);
//	stan::test::expect_ad(f_w(y, a, t0, v, sv, sw, st0), w);
//	stan::test::expect_ad(f_v(y, a, t0, w, sv, sw, st0), v);
//	stan::test::expect_ad(f_sv(y, a, t0, w, v, sw, st0), sv);
//	stan::test::expect_ad(f_sw(y, a, t0, w, v, sv, st0), sw);
//	stan::test::expect_ad(f_st0(y, a, t0, w, v, sv, sw), st0);



}



TEST(mathMixScalFun_w_sw_st0, wiener_full_lpdf) {
auto f_w_sw_st0 = [](const auto& y, const auto& a, const auto& t0, const auto& v, const auto& sv) {
     return [&y, &a, &t0, &v, &sv](const auto& w, const auto& sw, const auto& st0) {
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
      
      stan::test::expect_ad(f_w_sw_st0(y, a, t0, v, sv), w, sw, st0);
    }



TEST(mathMixScalFun_v_sv_sw, wiener_full_lpdf) {
auto f_v_sv_sw = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& st0) {
     return [&y, &a, &t0, &w, &st0](const auto& v, const auto& sv, const auto& sw) {
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
      
      stan::test::expect_ad(f_v_sv_sw(y, a, t0, w, st0), v, sv, sw);
    }



TEST(mathMixScalFun_v_sv_st0, wiener_full_lpdf) {
auto f_v_sv_st0 = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& sw) {
     return [&y, &a, &t0, &w, &sw](const auto& v, const auto& sv, const auto& st0) {
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
      
      stan::test::expect_ad(f_v_sv_st0(y, a, t0, w, sw), v, sv, st0);
    }



TEST(mathMixScalFun_v_sw_st0, wiener_full_lpdf) {
auto f_v_sw_st0 = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& sv) {
     return [&y, &a, &t0, &w, &sv](const auto& v, const auto& sw, const auto& st0) {
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
      
      stan::test::expect_ad(f_v_sw_st0(y, a, t0, w, sv), v, sw, st0);
    }



TEST(mathMixScalFun_sv_sw_st0, wiener_full_lpdf) {
auto f_sv_sw_st0 = [](const auto& y, const auto& a, const auto& t0, const auto& w, const auto& v) {
     return [&y, &a, &t0, &w, &v](const auto& sv, const auto& sw, const auto& st0) {
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
      
      stan::test::expect_ad(f_sv_sw_st0(y, a, t0, w, v), sv, sw, st0);
    }

