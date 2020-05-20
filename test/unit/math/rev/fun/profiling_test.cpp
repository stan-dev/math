#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <boost/random/mersenne_twister.hpp>

stan::math::profiles profiles;
TEST(Profiling, simple) {
    using stan::math::matrix_d;
    using stan::math::matrix_v;
    using stan::math::var;
    using stan::math::wishart_rng;
    using stan::math::cholesky_decompose;
    using stan::math::start_profiling;
    using stan::math::stop_profiling;

    boost::random::mt19937 rng;
    matrix_d I(800, 800);
    I.setZero();
    I.diagonal().setOnes();

    matrix_v Y = wishart_rng(5000, I, rng);
    matrix_v PP = Eigen::VectorXd::Ones(800);

    start_profiling(0, profiles, Y);
    matrix_v L = cholesky_decompose(Y);
    stop_profiling(0, profiles, L);

    start_profiling(1, profiles, L);
    matrix_v P = L * PP;
    stop_profiling(1, profiles, P);

    P(0,0).grad();

    std::chrono::duration<double> diff0f = profiles[0].fwd_pass_time_stop - profiles[0].fwd_pass_time_start;
    std::chrono::duration<double> diff0b = profiles[0].bkcwd_pass_time_stop - profiles[0].bkcwd_pass_time_start;
    std::chrono::duration<double> diff1f = profiles[1].fwd_pass_time_stop - profiles[1].fwd_pass_time_start;
    std::chrono::duration<double> diff1b = profiles[1].bkcwd_pass_time_stop - profiles[1].bkcwd_pass_time_start;

    std::cout << "0: FWD: " << diff0f.count() << ", BCKWD: " << diff0b.count() << std::endl;
    std::cout << "1: FWD: " << diff1f.count() << ", BCKWD: " << diff1b.count() << std::endl;
}
