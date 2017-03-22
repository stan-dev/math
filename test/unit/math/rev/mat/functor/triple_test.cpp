#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/triple.hpp>
#include <gtest/gtest.h>
#include <iostream>

TEST(MathMatrix, triple) {
  /* The jacobian is nested inside
   * the chain procedure.
   * NOTE: This test crashes at line 44 - 47 (third
   * iteration.)
   */
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;
  using stan::math::vari;
  using stan::math::ChainableStack;
  using stan::math::empty_nested;
  using stan::math::nested_size;
  using std::cout;
  using std::endl;

  int n_x = 3;

  for (int k = 0; k < n_x; k++) {
    Matrix<var, Dynamic, 1> x(n_x);
    x << 1, 2, 3;
    Matrix<var, Dynamic, 1> y = stan::math::triple(x);

    AVEC x_vec = createAVEC(x(0), x(1), x(2));
    VEC g;

    // Replicate grad function
      // Replicate void grad(vari* vi)
      typedef std::vector<vari*>::reverse_iterator it_t;
      y(k).vi_->init_dependent();
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

    MatrixXd Jf_x(n_x, n_x);
    Jf_x << 3, 0, 0, 0, 3, 0, 0, 0, 3;

    for (int i = 0; i < n_x; i++)
      EXPECT_EQ(Jf_x(k, i), g[i]);
  }
}

