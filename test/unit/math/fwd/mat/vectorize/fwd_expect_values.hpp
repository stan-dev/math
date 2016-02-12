#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_FWD_EXPECT_VALUES_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_FWD_EXPECT_VALUES_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <Eigen/Dense>

template <typename F>
static inline std::vector<double>
build_valid_fvar_vector(std::vector<double> double_vector,
                          int seed_index = -1) { 
  return F::valid_inputs();
}

template <typename F, typename V>
static inline std::vector<stan::math::fvar<V> >
build_valid_fvar_vector(std::vector<stan::math::fvar<V> > fvar_vector,
                          int seed_index = -1) { 
  using std::vector;
  using stan::math::fvar;

  vector<V> template_vector =
    build_valid_fvar_vector<F>(vector<V>(), seed_index);

  for (int i = 0; i < template_vector.size(); ++i) {
    if (seed_index == i)
      fvar_vector.push_back(
        fvar<V>(template_vector[i], template_vector[i]));
    else
      fvar_vector.push_back(fvar<V>(template_vector[i]));
  }
  return fvar_vector;
}

template <typename F, typename T>
static inline Eigen::Matrix<typename Eigen::internal::traits<T>::Scalar, 
                              T::RowsAtCompileTime,
                                T::ColsAtCompileTime>
build_valid_fvar_matrix(const T& x, int seed_index = -1) {

  typedef typename Eigen::internal::traits<T>::Scalar fvar_type;
  Eigen::Matrix<fvar_type, T::RowsAtCompileTime, T::ColsAtCompileTime>
    fvar_matrix(x.rows(), x.cols());
  size_t num_inputs = F::valid_inputs().size();
  for (int i = 0; i < x.size(); ++i) {
    std::vector<fvar_type> inputs;
    if (seed_index == i)
      inputs = build_valid_fvar_vector<F>(std::vector<fvar_type>(), 
                                       (seed_index % num_inputs));
    else
      inputs = build_valid_fvar_vector<F>(std::vector<fvar_type>()); 
    fvar_matrix(i) = inputs[(i % num_inputs)];
  }
  return fvar_matrix;
}

static inline void test_fvar(double exp_var, double test_var) {
  EXPECT_FLOAT_EQ(exp_var, test_var);
}

template <typename V>
static inline void test_fvar(V exp_var, V test_var) {
  test_fvar(exp_var.val(), test_var.val());
  test_fvar(exp_var.d_, test_var.d_);
}

template <typename F, typename V>
void expect_scalar_value() {
  using stan::math::fvar;
  using stan::test::expect_match_return_t;
  using std::vector;

  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<V> y = 
      build_valid_fvar_vector<F>(vector<V>(), i);
    V fy = F::template apply<V>(y[i]);
    V exp_y = F::apply_base(y[i]);
    test_fvar(exp_y, fy);
  }

  expect_match_return_t<V, V>();
  expect_match_return_t<std::vector<V>, 
                          std::vector<V> >();

}

template <typename F, typename V>
void expect_std_vector_values() {
  using stan::math::foo;
  using std::vector;
  using stan::math::fvar;
  using stan::math::foo;

  size_t num_inputs = F::valid_inputs().size();

  for (size_t i = 0; i < num_inputs; ++i) {
    vector<V> y =
      build_valid_fvar_vector<F>(vector<V>(), i);
    vector<V> fy = F::template apply<vector<V> >(y);
    EXPECT_EQ(y.size(), fy.size());
    test_fvar(F::apply_base(y[i]), fy[i]);
  }

  size_t vector_vector_size = 2;

  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {

      vector<vector<V> > z;
      for (size_t k = 0; k < num_inputs; ++k) {
        if (i == k)
          z.push_back(build_valid_fvar_vector<F>(vector<V>(), j));
        else
          z.push_back(build_valid_fvar_vector<F>(vector<V>()));
      }
      vector<vector<V> > fz = 
        F::template apply<vector<vector<V> > >(z);

      EXPECT_EQ(z.size(), fz.size());
      EXPECT_EQ(z[i].size(), fz[i].size());
      test_fvar(F::apply_base(z[i][j]), fz[i][j]);
    }
  }
}    

template <typename F, typename V>
void expect_matrix_values() {
  using stan::math::fvar;
  using std::vector;
  typedef 
    Eigen::Matrix<V, Eigen::Dynamic, Eigen::Dynamic> MatrixXvar;

  size_t num_inputs = F::valid_inputs().size();
  size_t num_cols = 3;
  MatrixXvar template_matrix(num_inputs, num_cols);

  for (int i = 0; i < template_matrix.size(); ++i) {
    MatrixXvar a = build_valid_fvar_matrix<F>(template_matrix, i);
    MatrixXvar fa = F::template apply<MatrixXvar>(a);
    EXPECT_EQ(a.size(), fa.size());
    test_fvar(F::apply_base(a(i)), fa(i));
  }

  size_t vector_matrix_size = 2;
  for (int i = 0; i < vector_matrix_size; ++i) {
    for (int j = 0; j < template_matrix.size(); ++j) {
      vector<MatrixXvar> b;
      for (int k = 0; k < vector_matrix_size; ++k)
        if (k == i)
          b.push_back(build_valid_fvar_matrix<F>(template_matrix, j));
        else
          b.push_back(build_valid_fvar_matrix<F>(template_matrix));
      vector<MatrixXvar> fb = F::template apply<vector<MatrixXvar> >(b);
      EXPECT_EQ(b.size(), fb.size());
      EXPECT_EQ(b[i].size(), fb[i].size());
      test_fvar(F::apply_base(b[i](j)), fb[i](j));
    }
  }

  int block_i = 1;
  int block_j = 1;
  int seed_i = block_j * num_inputs + block_i;
  MatrixXvar a = build_valid_fvar_matrix<F>(template_matrix, seed_i);
  MatrixXvar fab = foo(a.block(block_i, block_j, 1, 1));
  test_fvar(F::apply_base(a(1,1)), fab(0,0));
}

template <typename F, typename V>
void expect_vector_values() {
  using stan::math::fvar;
  using std::vector;
  typedef Eigen::Matrix<V, Eigen::Dynamic, 1> VectorXvar;

  size_t num_inputs = F::valid_inputs().size();
  VectorXvar template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    VectorXvar b = build_valid_fvar_matrix<F>(template_vector, i);
    VectorXvar fb = F::template apply<VectorXvar>(b);
    EXPECT_EQ(b.size(), fb.size());
    test_fvar(F::apply_base(b(i)), fb(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<VectorXvar> c;
      for (size_t k = 0; k < vector_vector_size; ++k)
        if (k == i)
          c.push_back(build_valid_fvar_matrix<F>(template_vector, j));
        else
          c.push_back(build_valid_fvar_matrix<F>(template_vector));
      vector<VectorXvar> fc = F::template apply<vector<VectorXvar> >(c);

      EXPECT_EQ(c.size(), fc.size());
      EXPECT_EQ(c[i].size(), fc[i].size());
      test_fvar(F::apply_base(c[i](j)), fc[i](j));
    }
  }
}

template <typename F, typename V>
void expect_row_vector_values() {
  using stan::math::fvar;
  using std::vector;
  typedef Eigen::Matrix<V, 1, Eigen::Dynamic> RowVectorXvar;

  size_t num_inputs = F::valid_inputs().size();
  RowVectorXvar template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    RowVectorXvar b = build_valid_fvar_matrix<F>(template_vector, i);
    RowVectorXvar fb = F::template apply<RowVectorXvar>(b);
    EXPECT_EQ(b.size(), fb.size());
    test_fvar(F::apply_base(b(i)), fb(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<RowVectorXvar> c;
      for (size_t k = 0; k < vector_vector_size; ++k)
        if (k == i) 
          c.push_back(build_valid_fvar_matrix<F>(template_vector, j));
        else
          c.push_back(build_valid_fvar_matrix<F>(template_vector));
      vector<RowVectorXvar> fc = 
          F::template apply<vector<RowVectorXvar> >(c);

      EXPECT_EQ(c.size(), fc.size());
      EXPECT_EQ(c[i].size(), fc[i].size());
      test_fvar(F::apply_base(c[i](j)), fc[i](j));
    }
  }
}

template <typename F>
void expect_values() {
  using stan::math::fvar;
  using std::vector;

  expect_scalar_value<F, fvar<double> >();
  expect_scalar_value<F, fvar<fvar<double> > >();
  expect_std_vector_values<F, fvar<double> >();
  expect_std_vector_values<F, fvar<fvar<double> > >();
  expect_matrix_values<F, fvar<double> >();
  expect_matrix_values<F, fvar<fvar<double> > >();
  expect_vector_values<F, fvar<double> >();
  expect_vector_values<F, fvar<fvar<double> > >();
  expect_row_vector_values<F, fvar<double> >();
  expect_row_vector_values<F, fvar<fvar<double> > >();
}

#endif
