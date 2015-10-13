#include <stan/math.hpp>
#include <stan/math/prim/scal/fun/selfstart_asymp.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, selfstart_asymp) {
  using stan::math::selfstart_asymp;

  double Asym = 100;
  double resp0 = -8.5;
  double lrc = -3.2;

  std::vector<double> input;
  input.push_back(0);
  input.push_back(-1);
  input.push_back(3.14);

  std::vector<double> expected_results;
  expected_results.push_back(-8.5);
  expected_results.push_back(-13.014075958479153883);
  expected_results.push_back(4.5352707075194018671);

  std::vector<double> test_results = 
    selfstart_asymp(input, Asym, resp0, lrc);

  EXPECT_EQ(expected_results.size(), test_results.size());
  for (std::size_t i = 0; i < test_results.size(); ++i) 
    EXPECT_FLOAT_EQ(expected_results[i], test_results[i]);  
}

TEST(MathFunctions, selfstart_asymp_error) {
  using stan::math::selfstart_asymp;

  double Asym = 100;
  double resp0 = -8.5;
  double lrc = -3.2;

  std::vector<double> input;
  input.push_back(0);
  input.push_back(-1);
  input.push_back(3.14);

  EXPECT_THROW(selfstart_asymp(input,stan::math::positive_infinity(), resp0,
                                 lrc),
               std::domain_error);
  EXPECT_THROW(selfstart_asymp(input, Asym, stan::math::positive_infinity(),
                                 lrc),
               std::domain_error);
  EXPECT_THROW(selfstart_asymp(input, Asym, resp0, 
                                 stan::math::positive_infinity()),
               std::domain_error);
  EXPECT_THROW(selfstart_asymp(input,stan::math::positive_infinity(), resp0,
                                 lrc),
               std::domain_error);
  input.push_back(std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW(selfstart_asymp(input, Asym, resp0, lrc), std::domain_error);
}

TEST(AgradRevMatrix, selfstart_asymp_rev) {
  using stan::math::selfstart_asymp;
  using stan::math::var;

  std::vector<double> input;
  input.push_back(3.14);
  
  var Asym = 100;
  double resp0 = -8.5;
  double lrc = -3.2;
  
  std::vector<var> f = selfstart_asymp(input, Asym, resp0, lrc);
  std::vector<var> params;
  params.push_back(Asym);  

  std::vector<double> x;
  f[0].grad(params,x);
  EXPECT_FLOAT_EQ(0.12014074384810513596, x[0]);
}
