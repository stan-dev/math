#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrixPower, two_by_two) {
  using stan::math::matrix_power;

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2, 2);
  Eigen::MatrixXd M(2, 2);
  M << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd M2(2, 2);
  M2 << 7.0, 10.0, 15.0, 22.0;
  Eigen::MatrixXd M3(2, 2);
  M3 << 37.0, 54.0, 81.0, 118.0;
  Eigen::MatrixXd M9(2, 2);
  M9 << 890461.0, 1297782.0, 1946673.0, 2837134.0;

  EXPECT_MATRIX_FLOAT_EQ(I, matrix_power(M, 0));
  EXPECT_MATRIX_FLOAT_EQ(M, matrix_power(M, 1));
  EXPECT_MATRIX_FLOAT_EQ(M2, matrix_power(M, 2));
  EXPECT_MATRIX_FLOAT_EQ(M3, matrix_power(M, 3));
  EXPECT_MATRIX_FLOAT_EQ(M9, matrix_power(M, 9));
}

TEST(MathMatrixPower, one_by_one) {
  using stan::math::matrix_power;

  Eigen::MatrixXd I(1, 1);
  I << 1.0;
  Eigen::MatrixXd M(1, 1);
  M << 3.0;
  Eigen::MatrixXd M2(1, 1);
  M2 << 9.0;
  Eigen::MatrixXd M3(1, 1);
  M3 << 27.0;
  Eigen::MatrixXd M9(1, 1);
  M9 << 19683.0;

  EXPECT_MATRIX_FLOAT_EQ(I, matrix_power(M, 0));
  EXPECT_MATRIX_FLOAT_EQ(M, matrix_power(M, 1));
  EXPECT_MATRIX_FLOAT_EQ(M2, matrix_power(M, 2));
  EXPECT_MATRIX_FLOAT_EQ(M3, matrix_power(M, 3));
  EXPECT_MATRIX_FLOAT_EQ(M9, matrix_power(M, 9));
}

TEST(MathMatrixPower, compare_to_simple_impl) {
  using stan::math::matrix_power;

  int size = 5;
  Eigen::MatrixXd M = Eigen::MatrixXd::Random(size, size);
  Eigen::MatrixXd expected = Eigen::MatrixXd::Identity(size, size);
  int exponent = 4;
  for (int i = 0; i < exponent; i++)
    expected *= M;

  EXPECT_MATRIX_FLOAT_EQ(expected, matrix_power(M, exponent));
}

TEST(MathMatrixPower, zero_size) {
  using stan::math::matrix_power;
  Eigen::MatrixXd zero_size = Eigen::MatrixXd::Identity(0, 0);
  EXPECT_NO_THROW(matrix_power(zero_size, 2));
}

TEST(MathMatrixPower, not_square) {
  using stan::math::matrix_power;
  Eigen::MatrixXd not_square = Eigen::MatrixXd::Identity(3, 4);
  EXPECT_THROW(matrix_power(not_square, 2), std::invalid_argument);
}

TEST(MathMatrixPower, negative_exponent) {
  using stan::math::matrix_power;
  Eigen::MatrixXd good = Eigen::MatrixXd::Identity(2, 2);
  EXPECT_NO_THROW(matrix_power(good, 2));
  EXPECT_THROW(matrix_power(good, -2), std::domain_error);
}

TEST(MathMatrixPower, inf_and_nan) {
  using stan::math::matrix_power;

  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  double ninf = -inf;

  Eigen::MatrixXd good = Eigen::MatrixXd::Identity(2, 2);
  EXPECT_NO_THROW(matrix_power(good, 2));

  good(0, 0) = nan;
  EXPECT_THROW(matrix_power(good, 2), std::domain_error);
  good(0, 0) = inf;
  EXPECT_THROW(matrix_power(good, 2), std::domain_error);
  good(0, 0) = ninf;
  EXPECT_THROW(matrix_power(good, 2), std::domain_error);
}

TEST(MathMatrixPower, invalid_argument_vs_domain_error) {
  using stan::math::matrix_power;
  Eigen::MatrixXd not_square = Eigen::MatrixXd::Identity(3, 4);
  Eigen::MatrixXd good = Eigen::MatrixXd::Identity(2, 2);
  double nan = std::numeric_limits<double>::quiet_NaN();
  not_square(0, 0) = nan;
  good(0, 0) = nan;
  EXPECT_THROW(matrix_power(not_square, 2), std::invalid_argument);
  EXPECT_THROW(matrix_power(good, -2), std::domain_error);
}

TEST(MathMatrixPower, matrix_power_operator) {
  using stan::math::operator^;

  Eigen::MatrixXd M(2, 2);
  M << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd M2(2, 2);
  M2 << 7.0, 10.0, 15.0, 22.0;

  EXPECT_MATRIX_FLOAT_EQ(M2, (M ^ 2));
}

TEST(MathMatrixPower, large) {
  using stan::math::matrix_power;

  Eigen::MatrixXd M(10, 10);
  M << 0.85704752, 0.00243432, 0.56797957, 0.09797887, 0.79268669, 0.33956719,
      0.13723404, 0.84215282, 0.28619714, 0.08314874, 0.07948192, 0.94349737,
      0.29626239, 0.66144137, 0.80467646, 0.66725623, 0.75609966, 0.71926087,
      0.13438711, 0.20910412, 0.09590056, 0.14401656, 0.88054535, 0.9165302,
      0.23232802, 0.6012156, 0.20750787, 0.16975039, 0.09254274, 0.48880858,
      0.97982399, 0.76057945, 0.51750248, 0.80098624, 0.5320405, 0.81460258,
      0.65338107, 0.59966065, 0.80510786, 0.58025303, 0.39742551, 0.60716168,
      0.29986078, 0.69634038, 0.3687128, 0.35146949, 0.82494471, 0.5249718,
      0.86960285, 0.23226187, 0.5070069, 0.21085213, 0.24510391, 0.4985713,
      0.69277847, 0.73104199, 0.40068102, 0.28178094, 0.36478028, 0.76132117,
      0.68874095, 0.25514365, 0.06871063, 0.74883898, 0.41938939, 0.426707,
      0.5631822, 0.54574921, 0.36608806, 0.18924788, 0.27345123, 0.29301448,
      0.51680293, 0.37062514, 0.30082761, 0.70411961, 0.45501945, 0.08786704,
      0.82768645, 0.30843809, 0.27834481, 0.43559884, 0.16661187, 0.82504191,
      0.04386024, 0.96918174, 0.97039513, 0.30973888, 0.50857728, 0.89804354,
      0.59800669, 0.85096666, 0.02926879, 0.41940341, 0.65511543, 0.16498608,
      0.44015271, 0.34145677, 0.56585428, 0.46709031;

  Eigen::MatrixXd M9(10, 10);
  M9 << 126466.445764333, 112623.642622401, 86368.2907351496, 150340.210416751,
      123356.489363506, 146624.574135454, 137715.026747634, 113774.387130492,
      124140.892158493, 107237.515347083, 187371.757177332, 166862.308458682,
      127962.720718795, 222743.051924318, 182764.054690347, 217238.052371467,
      204037.749772744, 168567.371596176, 183926.388050216, 158882.338085627,
      137288.259112996, 122260.997166868, 93759.0543967740, 163205.152878243,
      133912.244199050, 159171.600155708, 149499.689686175, 123510.288692996,
      134763.926331194, 116413.944342069, 242832.305977762, 216252.299976024,
      165838.772088294, 288673.424986375, 236860.779441396, 281538.959245798,
      264431.432537886, 218462.001221056, 238367.205906074, 205910.345107489,
      183522.327745970, 163434.291897440, 125333.910074745, 218167.023591516,
      179009.341446653, 212775.121253318, 199846.033154285, 165104.348258986,
      180147.828488179, 155618.207131273, 162117.333177683, 144372.231603063,
      110715.677560648, 192721.338773521, 158130.675391508, 187958.307184737,
      176537.149432354, 145847.477005024, 159136.408052580, 137467.855177878,
      149780.220252025, 133385.538453015, 102290.221233416, 178055.276643238,
      146096.944642611, 173654.685223753, 163102.653665707, 134748.470157240,
      147026.106115242, 127006.588278837, 144142.315282588, 128364.773959074,
      98439.9365839541, 171353.006608093, 140597.754128395, 167118.106131280,
      156963.357862540, 129676.501068110, 141491.963957078, 122225.804438213,
      190602.006087896, 169739.052954256, 130168.979099968, 226583.382876596,
      185914.907329016, 220983.481170086, 207555.590947477, 171473.548964940,
      187097.435673739, 161621.544395111, 160013.515879562, 142498.647686615,
      109278.842978772, 190220.314850057, 156078.552041347, 185519.091744448,
      174246.126849153, 143954.733727807, 157071.213262326, 135683.920871181;

  EXPECT_MATRIX_FLOAT_EQ(M9, matrix_power(M, 9));
}
