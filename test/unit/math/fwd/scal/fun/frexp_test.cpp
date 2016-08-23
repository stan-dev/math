#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdFrexp, Fvar) {
    using stan::math::fvar;
    using std::frexp;
    
    fvar<double> x(123.45, 1.0);
    fvar<double> y(8.0, 2.0);
    int alpha, beta, gamma, delta, epsilon, zeta;
    
    fvar<double> a = frexp(x, &alpha);
    EXPECT_FLOAT_EQ(frexp(123.45, &beta), a.val_);
    EXPECT_EQ(alpha, beta);
    EXPECT_FLOAT_EQ(0.0, a.d_);

    fvar<double> b = frexp(y, &gamma);
    EXPECT_FLOAT_EQ(frexp(8.0, &delta), b.val_);
    EXPECT_EQ(gamma, delta);
    EXPECT_FLOAT_EQ(0.0, b.d_);
    
    fvar<double> c = frexp(2 * x, &epsilon);
    EXPECT_FLOAT_EQ(frexp(2 * 123.45, &zeta), c.val_);
    EXPECT_EQ(epsilon, zeta);
    EXPECT_FLOAT_EQ(0.0, c.d_);
}

TEST(AGradFwdFrexp, FvarFvarDouble) {
    using stan::math::fvar;
    using std::frexp;
    
    int alpha, beta;
    
    fvar<fvar<double> > x;
    x.val_.val_ = 16.0;
    x.val_.d_ = 2.0;
    
    fvar<fvar<double> > a = frexp(x, &alpha);
    
    EXPECT_FLOAT_EQ(frexp(16.0, &beta), a.val_.val_);
    EXPECT_EQ(alpha, beta);
    EXPECT_FLOAT_EQ(0, a.val_.d_);
    EXPECT_FLOAT_EQ(0, a.d_.val_);
    EXPECT_FLOAT_EQ(0, a.d_.d_);
    
    
    fvar<fvar<double> > y;
    y.val_.val_ = 16.0;
    y.d_.val_ = 2.0;
    
  //  a = frexp(y, &alpha);
    EXPECT_FLOAT_EQ(frexp(16.0, &beta), a.val_.val_);
    EXPECT_FLOAT_EQ(0, a.val_.d_);
    EXPECT_FLOAT_EQ(0, a.d_.val_);
    EXPECT_FLOAT_EQ(0, a.d_.d_);
}

struct frexp_fun {
    template <typename T0>
    inline T0
    operator()(const T0& arg1) const {
        int alpha;
        return frexp(arg1, &alpha);
    }
};

TEST(AgradFwdFrexp, frexp_NaN) {
    frexp_fun frexp_;
    test_nan_fwd(frexp_,false); // CHECK - no exception (second argument)?
}
