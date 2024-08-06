#ifndef STAN_TEST_UNIT_MIX_LAPLACE_MOTORCYCLE_GP_HPP
#define STAN_TEST_UNIT_MIX_LAPLACE_MOTORCYCLE_GP_HPP
namespace stan {
namespace test {
namespace laplace {
namespace moto {
auto x = std::vector<double>{
    2.4,  2.6,  3.2,  3.6,  4,    6.2,  6.6,  6.8,  7.8,  8.2,  8.8,  8.8,
    9.6,  10,   10.2, 10.6, 11,   11.4, 13.2, 13.6, 13.8, 14.6, 14.6, 14.6,
    14.6, 14.6, 14.6, 14.8, 15.4, 15.4, 15.4, 15.4, 15.6, 15.6, 15.8, 15.8,
    16,   16,   16.2, 16.2, 16.2, 16.4, 16.4, 16.6, 16.8, 16.8, 16.8, 17.6,
    17.6, 17.6, 17.6, 17.8, 17.8, 18.6, 18.6, 19.2, 19.4, 19.4, 19.6, 20.2,
    20.4, 21.2, 21.4, 21.8, 22,   23.2, 23.4, 24,   24.2, 24.2, 24.6, 25,
    25,   25.4, 25.4, 25.6, 26,   26.2, 26.2, 26.4, 27,   27.2, 27.2, 27.2,
    27.6, 28.2, 28.4, 28.4, 28.6, 29.4, 30.2, 31,   31.2, 32,   32,   32.8,
    33.4, 33.8, 34.4, 34.8, 35.2, 35.2, 35.4, 35.6, 35.6, 36.2, 36.2, 38,
    38,   39.2, 39.4, 40,   40.4, 41.6, 41.6, 42.4, 42.8, 42.8, 43,   44,
    44.4, 45,   46.6, 47.8, 47.8, 48.8, 50.6, 52,   53.2, 55,   55,   55.4,
    57.6};
auto y = Eigen::VectorXd{
    {0,      -1.3,   -2.7,   0,      -2.7,   -2.7,   -2.7,   -1.3,   -2.7,
     -2.7,   -1.3,   -2.7,   -2.7,   -2.7,   -5.4,   -2.7,   -5.4,   0,
     -2.7,   -2.7,   0,      -13.3,  -5.4,   -5.4,   -9.3,   -16,    -22.8,
     -2.7,   -22.8,  -32.1,  -53.5,  -54.9,  -40.2,  -21.5,  -21.5,  -50.8,
     -42.9,  -26.8,  -21.5,  -50.8,  -61.7,  -5.4,   -80.4,  -59,    -71,
     -91.1,  -77.7,  -37.5,  -85.6,  -123.1, -101.9, -99.1,  -104.4, -112.5,
     -50.8,  -123.1, -85.6,  -72.3,  -127.2, -123.1, -117.9, -134,   -101.9,
     -108.4, -123.1, -123.1, -128.5, -112.5, -95.1,  -81.8,  -53.5,  -64.4,
     -57.6,  -72.3,  -44.3,  -26.8,  -5.4,   -107.1, -21.5,  -65.6,  -16,
     -45.6,  -24.2,  9.5,    4,      12,     -21.5,  37.5,   46.9,   -17.4,
     36.2,   75,     8.1,    54.9,   48.2,   46.9,   16,     45.6,   1.3,
     75,     -16,    -54.9,  69.6,   34.8,   32.1,   -37.5,  22.8,   46.9,
     10.7,   5.4,    -1.3,   -21.5,  -13.3,  30.8,   -10.7,  29.4,   0,
     -10.7,  14.7,   -1.3,   0,      10.7,   10.7,   -26.8,  -14.7,  -13.3,
     0,      10.7,   -14.7,  -2.7,   10.7,   -2.7,   10.7}};
}  // namespace moto
}  // namespace laplace
}  // namespace test
}  // namespace stan
#endif
