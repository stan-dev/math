#ifndef TEST_UNIT_MATH_PRIM_FUNCTOR_ODE_FIXTURE
#define TEST_UNIT_MATH_PRIM_FUNCTOR_ODE_FIXTURE

struct StanMathOde : public ::testing::Test {
  std::stringstream msgs;
  std::vector<double> x;
  std::vector<int> x_int;
};

#endif
