#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/triple.hpp>
#include <gtest/gtest.h>
#include <iostream>

template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1>
triple_f(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) {
  int nStates = 3;
  Eigen::Matrix<T, Eigen::Dynamic, 1> y(nStates);
  y(0) = 3 * x(0);
  y(1) = 3 * x(1);
  y(2) = 3 * x(2);
  return y;
}

struct triple_f_functor {
public:
  triple_f_functor() { };

  template <typename T>
  inline
  Eigen::Matrix<T, Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
    return triple_f(x);
  }
};

TEST(MathMatrix, jacobian_1) {
  /* TEST 1: the jacobian is computed after
   * we compute y = f(x_var) but before we run 
   * grad.
   * Note: This test works.
   */
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;

  int n_x = 3;
  for (int k = 0; k < n_x; k++) {
    VectorXd x(n_x);
    x << 1, 2, 3;
    Matrix<var, Dynamic, 1> x_var(n_x);
    for (int i = 0; i < n_x; i++)
      x_var(i) = x(i);

    triple_f_functor f;
    Matrix<var, Dynamic, 1> y = f(x_var);

    VectorXd fvec;
    MatrixXd J;
    stan::math::jacobian(f, stan::math::value_of(x), fvec, J);

    MatrixXd Jf_x(n_x, n_x);
    Jf_x << 3, 0, 0, 0, 3, 0, 0, 0, 3;

    for (int i = 0; i < n_x; i++)
      for (int j = 0; j < n_x; j++)
        EXPECT_EQ(Jf_x(i, j), J(i, j));

    VectorXd solutions(n_x);
    solutions << 3, 6, 9;

    EXPECT_EQ(solutions(k), y(k));

    AVEC x_vec = createAVEC(x_var(0), x_var(1), x_var(2));
    VEC g;
    y(k).grad(x_vec, g);

    for (int i = 0; i < n_x; i++)
      EXPECT_EQ(Jf_x(k, i), g[i]);
  }
}  

TEST(MathMatrix, jacobian_2) {
  /* TEST 2: the jacobian is nested inside
   * the grad procedure.
   * Note: This test works.
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

  int n_x = 3;

  for (int k = 0; k < n_x; k++) {
    VectorXd x(n_x);
    x << 1, 2, 3;
    Matrix<var, Dynamic, 1> x_var(n_x);
    for (int i = 0; i < n_x; i++)
      x_var(i) = x(i);

    triple_f_functor f;
    Matrix<var, Dynamic, 1> y = f(x_var);

    VectorXd solutions(n_x);
    solutions << 3, 6, 9;

    EXPECT_EQ(solutions(k), y(k));

    VectorXd fvec;
    MatrixXd J;

    AVEC x_vec = createAVEC(x_var(0), x_var(1), x_var(2));
    VEC g;

    // Replicate grad function
      // Replicate void grad(vari* vi)
      typedef std::vector<vari*>::reverse_iterator it_t;
      y(k).vi_->init_dependent();
      it_t begin = ChainableStack::var_stack_.rbegin();
      it_t end = empty_nested()
        ? ChainableStack::var_stack_.rend() : begin + nested_size();
      for (it_t it = begin; it < end; ++it) {
        (*it)->chain();
        stan::math::jacobian(f, stan::math::value_of(x), fvec, J);
      }

    g.resize(x_vec.size());
    for (size_t i = 0; i < x_vec.size(); ++i)
      g[i] = x_vec[i].vi_->adj_;

    MatrixXd Jf_x(n_x, n_x);
    Jf_x << 3, 0, 0, 0, 3, 0, 0, 0, 3;

    for (int i = 0; i < n_x; i++)
      EXPECT_EQ(Jf_x(k, i), g[i]);
  }
} 


TEST(MathMatrix, jacobian_3) {
  /* TEST 3: the jacobian is nested inside
   * the chain procedure.
   * NOTE: This test crashes.
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
        // cout << "o " << endl;  // shows for loop runs
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

