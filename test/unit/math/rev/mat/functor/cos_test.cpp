#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/cos.hpp>
#include <gtest/gtest.h>
#include <iostream>

TEST(JacobianTest, cos) {
  using stan::math::var;
  using stan::math::vari;
  using stan::math::ChainableStack;
  using stan::math::empty_nested;
  using stan::math::nested_size;
  using std::cout;
  using std::endl;

  var x = 0;
  var y = stan::math::cos(x);
  
  AVEC x_vec = createAVEC(x);
  VEC g;
  
  // Replicate grad function
    // Replicate void(vari* vi)
    typedef std::vector<vari*>::reverse_iterator it_t;
      y.vi_->init_dependent();
      it_t begin = ChainableStack::var_stack_.rbegin();
      it_t end = empty_nested()
        ? ChainableStack::var_stack_.rend() : begin + nested_size();
      cout << "Marker A" << endl;      
      for (it_t it = begin; it < end; ++it) {
        cout << "o " << endl;  // shows for loop runs
        (*it)->chain();  // the chain() has Jacobian nested inside
      }
      std::cout << "Marker B" << std::endl;
    
    g.resize(x_vec.size());
    for (size_t i = 0; i < x_vec.size(); ++i)
      g[i] = x_vec[i].vi_->adj_;
    
  EXPECT_EQ(1, y);
  EXPECT_EQ(0, g[0]);
}
