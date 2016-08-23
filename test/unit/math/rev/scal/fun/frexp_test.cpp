#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>

TEST(AgradRev, frexp) {
    AVAR a = 123.45;
    int b;
    AVAR f = frexp(a, &b);
    EXPECT_FLOAT_EQ(0.964453, f.val());
    EXPECT_EQ(7, b);
    
    AVEC x = createAVEC(a);
    VEC grad_f;
    f.grad(x, grad_f);
    EXPECT_FLOAT_EQ(0.0, grad_f[0]);
}

TEST(AgradRev, frexp_2){
    AVAR a = -123.45;
    int b;
    AVAR f = frexp(a, &b);
    EXPECT_FLOAT_EQ(-0.964453, f.val());
    EXPECT_EQ(7, b);
    
    AVEC x = createAVEC(a);
    VEC grad_f;
    f.grad(x, grad_f);
    EXPECT_FLOAT_EQ(0.0, grad_f[0]);
}

struct frexp_fun {
    template <typename T0>
    inline T0
    operator()(const T0& arg1) const {
        int b;
        return frexp(arg1, &b);
    }
};

TEST(AgradRev, frexp_NAN) {
    frexp_fun frexp_;
    test_nan(frexp_,false,true);
}