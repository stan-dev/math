#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_VALUES_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_VALUES_HPP

#include <stan/math/rev/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_match_return_t.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

template <typename F, typename T>
static inline Eigen::Matrix<stan::math::var, T::RowsAtCompileTime,
                              T::ColsAtCompileTime>
build_valid_var_matrix(const T& x) {
  Eigen::Matrix<stan::math::var, T::RowsAtCompileTime, 
    T::ColsAtCompileTime> var_matrix(x.rows(), x.cols());
  std::vector<double> inputs = F::valid_inputs();
  for (int i = 0; i < x.size(); ++i) {
      var_matrix(i) = inputs[(i % inputs.size())];
  }
  return var_matrix;
}

template <typename F>
static inline std::vector<stan::math::var>
build_valid_var_vector() {

  std::vector<double> inputs = F::valid_inputs();
  std::vector<stan::math::var> var_vector;

  for (int i = 0; i < inputs.size(); ++i) {
    var_vector.push_back(inputs[i]);
  }
  
  return var_vector;
}

template <typename F>
static inline void test_autodiff(double test_val, double test_adj) {
  using stan::math::var;

  var x = test_val;
  F::apply_base(x).grad();
  EXPECT_FLOAT_EQ(x.adj(), test_adj);
}

template <typename F>
void expect_scalar_value() {
  using stan::math::var;

  std::vector<var> y = build_valid_var_vector<F>();
  for (size_t i = 0; i < y.size(); ++i) {
    var fy = F::template apply<var>(y[i]);

    EXPECT_FLOAT_EQ(F::apply_base(y[i]).val(), fy.val());

    fy.grad();
    test_autodiff<F>(y[i].val(), y[i].adj());
  }
}

template <typename F>
void expect_std_vectors_value() {
  using std::vector;
  using stan::math::var;

  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    
    vector<var> y = build_valid_var_vector<F>();
    vector<var> fy = F::template apply<vector<var> >(y);

    EXPECT_EQ(y.size(), fy.size());
    EXPECT_FLOAT_EQ(F::apply_base(y[i]).val(), fy[i].val());

    fy[i].grad();
    test_autodiff<F>(y[i].val(), y[i].adj());
  }
}

template <typename F>
void expect_std_vector_vectors_value() {
  using std::vector;
  using stan::math::var;

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {

    for (size_t j = 0; j < F::valid_inputs().size(); ++j) {

      vector<vector<var> > z;
      for (size_t i = 0; i < vector_vector_size; ++i) {
        z.push_back(build_valid_var_vector<F>());
      } 
      vector<vector<var> > fz = 
                        F::template apply<vector<vector<var> > >(z);

      EXPECT_EQ(z.size(), fz.size());
      EXPECT_EQ(z[i].size(), fz[i].size());
      EXPECT_FLOAT_EQ(F::apply_base(z[i][j]).val(), fz[i][j].val());

      fz[i][j].grad();
      test_autodiff<F>(z[i][j].val(), z[i][j].adj());
    }
  }
}    

template <typename F>
void expect_matrix_value() {
  using stan::math::var;
  using std::vector;
  typedef Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> MatrixXvar;
  
  size_t num_cols = 3;
  size_t num_inputs = F::valid_inputs().size();
  MatrixXvar template_matrix(num_inputs, num_cols);

  for (size_t i = 0; i < template_matrix.size(); ++i) {
    MatrixXvar a = build_valid_var_matrix<F>(template_matrix);
    MatrixXvar fa = F::template apply<MatrixXvar>(a);
    EXPECT_EQ(a.size(), fa.size());
    EXPECT_FLOAT_EQ(F::apply_base(a(i)).val(), fa(i).val());

    fa(i).grad();
    test_autodiff<F>(a(i).val(), a(i).adj());
  }

  size_t vector_matrix_size = 2;
  for (size_t i = 0; i < vector_matrix_size; ++i) {
    for (size_t j = 0; j < template_matrix.size(); ++j) {

      vector<MatrixXvar> b;
      for (size_t k = 0; k < vector_matrix_size; ++k)
        b.push_back(build_valid_var_matrix<F>(template_matrix));
      vector<MatrixXvar> fb = F::template apply<vector<MatrixXvar> >(b);

      EXPECT_EQ(b[i].size(), fb[i].size());
      EXPECT_EQ(b[i].rows(), fb[i].rows());
      EXPECT_EQ(b[i].cols(), fb[i].cols());
      EXPECT_FLOAT_EQ(F::apply_base(b[i](j)).val(), fb[i](j).val()); 

      fb[i](j).grad();
      test_autodiff<F>(b[i](j).val(), b[i](j).adj());
    }
  }

  MatrixXvar a = build_valid_var_matrix<F>(template_matrix);
  MatrixXvar fab = F::template apply<MatrixXvar>(a.block(1, 1, 1, 1));
  EXPECT_FLOAT_EQ(F::apply_base(a(1,1)).val(), fab(0,0).val());
  fab(0,0).grad();
  test_autodiff<F>(a(1,1).val(), a(1,1).adj());
}

template <typename F>
void expect_vector_value() {
  using stan::math::var;
  using std::vector;
  typedef Eigen::Matrix<var, Eigen::Dynamic, 1> VectorXvar;

  size_t num_inputs = F::valid_inputs().size();
  VectorXvar template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    VectorXvar b = build_valid_var_matrix<F>(template_vector);
    VectorXvar fb = F::template apply<VectorXvar>(b);
    EXPECT_EQ(b.size(), fb.size());
    
    EXPECT_FLOAT_EQ(F::apply_base(b(i)).val(), fb(i).val());

    fb(i).grad();
    test_autodiff<F>(b(i).val(), b(i).adj());
  }

  size_t vector_vector_size = 2;

  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<VectorXvar> vb;
      for (size_t k = 0; k < num_inputs; ++k)
        vb.push_back(build_valid_var_matrix<F>(template_vector));
      vector<VectorXvar> fvb = F::template apply<vector<VectorXvar> >(vb);
      EXPECT_EQ(vb[i].size(), fvb[i].size());
      EXPECT_EQ(vb.size(), fvb.size());
      EXPECT_FLOAT_EQ(F::apply_base(vb[i](j)).val(), fvb[i](j).val());

      fvb[i](j).grad();
      test_autodiff<F>(vb[i](j).val(), vb[i](j).adj()); 
    }
  }
}

template <typename F>
void expect_row_vector_value() {
  using stan::math::var;
  using std::vector;
  typedef Eigen::Matrix<var, 1, Eigen::Dynamic> RowVectorXvar;

  size_t num_inputs = F::valid_inputs().size();
  RowVectorXvar template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    RowVectorXvar b = build_valid_var_matrix<F>(template_vector);
    RowVectorXvar fb = F::template apply<RowVectorXvar>(b);
    EXPECT_EQ(b.size(), fb.size());
    
    EXPECT_FLOAT_EQ(F::apply_base(b(i)).val(), fb(i).val());

    fb(i).grad();
    test_autodiff<F>(b(i).val(), b(i).adj());
  }

  size_t vector_vector_size = 2;

  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<RowVectorXvar> vb;
      for (size_t k = 0; k < num_inputs; ++k)
        vb.push_back(build_valid_var_matrix<F>(template_vector));
      vector<RowVectorXvar> fvb = 
        F::template apply<vector<RowVectorXvar> >(vb);
      EXPECT_EQ(vb[i].size(), fvb[i].size());
      EXPECT_EQ(vb.size(), fvb.size());
      EXPECT_FLOAT_EQ(F::apply_base(vb[i](j)).val(), fvb[i](j).val());

      fvb[i](j).grad();
      test_autodiff<F>(vb[i](j).val(), vb[i](j).adj()); 
    }
  }
}

template <typename F>
void expect_values() {
  expect_scalar_value<F>();
  expect_std_vectors_value<F>();
  expect_std_vector_vectors_value<F>();
  expect_matrix_value<F>();
  expect_vector_value<F>();
  expect_row_vector_value<F>();
}

#endif
