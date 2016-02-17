#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_MIX_EXPECT_VALUES_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_MIX_EXPECT_VALUES_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/fwd/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <Eigen/Dense>
#include <test/unit/math/rev/mat/fun/util.hpp>

template <typename F>
static inline std::vector<stan::math::var>
build_valid_fvar_vector(std::vector<stan::math::var> var_vector,
                            int seed_index = -1) { 
  std::vector<double> inputs = F::valid_inputs();
  for (size_t i = 0; i < inputs.size(); ++i) {
      var_vector.push_back(inputs[i]);
  }
  return var_vector;
}

template <typename F, typename V>
static inline std::vector<stan::math::fvar<V> >
build_valid_fvar_vector(std::vector<stan::math::fvar<V> > fvar_vector,
                          int seed_index = -1) { 
  using std::vector;
  using stan::math::fvar;

  vector<V> val_vector =
    build_valid_fvar_vector<F>(vector<V>(), seed_index);
  vector<V> d_vector; 
  if (seed_index != -1)
    d_vector = build_valid_fvar_vector<F>(vector<V>(), seed_index);

  for (size_t i = 0; i < val_vector.size(); ++i) {
    if (seed_index == static_cast<int>(i))
      fvar_vector.push_back(fvar<V>(val_vector[i], d_vector[i]));
    else
      fvar_vector.push_back(fvar<V>(val_vector[i]));
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

static inline void test_fvar(stan::math::var exp_var,
                               stan::math::var base_exp_var,
                                 stan::math::var test_var,
                                   stan::math::var base_test_var) {
  EXPECT_FLOAT_EQ(exp_var.val(), test_var.val());
  AVEC exp_y = createAVEC(base_exp_var);
  VEC exp_g;
  exp_var.grad(exp_y, exp_g);
  AVEC test_y = createAVEC(base_test_var);
  VEC test_g;
  test_var.grad(test_y, test_g);
  EXPECT_FLOAT_EQ(exp_g[0], test_g[0]);
}

template <typename V>
static inline void test_fvar(V exp_var, V base_exp_var, 
                               V test_var, V base_test_var) {
  test_fvar(exp_var.val(), base_exp_var.val(),
              test_var.val(), base_test_var.val());
  test_fvar(exp_var.d_, base_exp_var.d_, 
              test_var.d_, base_test_var.d_);
}

template <typename F, typename V>
void expect_scalar_value() {
  using stan::math::fvar;
  using stan::test::expect_match_return_t;
  using std::vector;
  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<V> y = build_valid_fvar_vector<F>(vector<V>(), i);
    vector<V> z = build_valid_fvar_vector<F>(vector<V>(), i);
    V fz = F::template apply<V>(z[i]);
    test_fvar(F::apply_base(y[i]), y[i], fz, z[i]);
  }
  expect_match_return_t<F, V, V>();
  expect_match_return_t<F, std::vector<V>, std::vector<V> >();
}

template <typename F, typename V>
void expect_std_vector_values() {
  using std::vector;
  using stan::math::fvar;

  size_t num_inputs = F::valid_inputs().size();

  for (size_t i = 0; i < num_inputs; ++i) {
    vector<V> y = build_valid_fvar_vector<F>(vector<V>(), i);
    vector<V> z = build_valid_fvar_vector<F>(vector<V>(), i);
    vector<V> fz = F::template apply<vector<V> >(z);
    EXPECT_EQ(z.size(), fz.size());
    test_fvar(F::apply_base(y[i]), y[i], fz[i], z[i]);
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<vector<V> > a;
      vector<vector<V> > b;
      for (size_t k = 0; k < num_inputs; ++k) {
        if (i == k) {
          a.push_back(build_valid_fvar_vector<F>(vector<V>(), j));
          b.push_back(build_valid_fvar_vector<F>(vector<V>(), j));
        }
        else {
          a.push_back(build_valid_fvar_vector<F>(vector<V>()));
          b.push_back(build_valid_fvar_vector<F>(vector<V>()));
        }
      }
      vector<vector<V> > fb = 
        F::template apply<vector<vector<V> > >(b);

      EXPECT_EQ(b.size(), fb.size());
      EXPECT_EQ(b[i].size(), fb[i].size());
      test_fvar(F::apply_base(a[i][j]), a[i][j], fb[i][j], b[i][j]);
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
    MatrixXvar y = build_valid_fvar_matrix<F>(template_matrix, i);
    MatrixXvar z = build_valid_fvar_matrix<F>(template_matrix, i);
    MatrixXvar fz = F::template apply<MatrixXvar>(z);
    EXPECT_EQ(z.size(), fz.size());
    test_fvar(F::apply_base(y(i)), y(i), fz(i), z(i));
  }

  size_t vector_matrix_size = 2;
  for (size_t i = 0; i < vector_matrix_size; ++i) {
    for (int j = 0; j < template_matrix.size(); ++j) {
      vector<MatrixXvar> a;
      vector<MatrixXvar> b;
      for (size_t k = 0; k < vector_matrix_size; ++k) {
        if (k == i) {
          a.push_back(build_valid_fvar_matrix<F>(template_matrix, j));
          b.push_back(build_valid_fvar_matrix<F>(template_matrix, j));
        } else {
          a.push_back(build_valid_fvar_matrix<F>(template_matrix));
          b.push_back(build_valid_fvar_matrix<F>(template_matrix));
        }
      }
      vector<MatrixXvar> fb = F::template apply<vector<MatrixXvar> >(b);
      EXPECT_EQ(b.size(), fb.size());
      EXPECT_EQ(b[i].size(), fb[i].size());
      test_fvar(F::apply_base(a[i](j)), a[i](j), fb[i](j), b[i](j));
    }
  }

  int block_i = 1;
  int block_j = 1;
  int seed_i = block_j * num_inputs + block_i; 
  MatrixXvar c = build_valid_fvar_matrix<F>(template_matrix, seed_i);
  MatrixXvar d = build_valid_fvar_matrix<F>(template_matrix, seed_i);
  MatrixXvar fab = foo(d.block(block_i, block_j, 1, 1));
  test_fvar(F::apply_base(c(block_i, block_j)), c(block_i, block_j),
              fab(0,0), d(block_i, block_j));
}

template <typename F, typename V>
void expect_vector_values() {
  using stan::math::fvar;
  using std::vector;
  typedef Eigen::Matrix<V, Eigen::Dynamic, 1> VectorXvar;

  size_t num_inputs = F::valid_inputs().size();
  VectorXvar template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    VectorXvar a = build_valid_fvar_matrix<F>(template_vector, i);
    VectorXvar b = build_valid_fvar_matrix<F>(template_vector, i);
    VectorXvar fb = F::template apply<VectorXvar>(b);
    EXPECT_EQ(b.size(), fb.size());
    test_fvar(F::apply_base(a(i)), a(i), fb(i), b(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<VectorXvar> c;
      vector<VectorXvar> d;
      for (size_t k = 0; k < vector_vector_size; ++k)
        if (k == i) {
          c.push_back(build_valid_fvar_matrix<F>(template_vector, j));
          d.push_back(build_valid_fvar_matrix<F>(template_vector, j));
        }
        else {
          c.push_back(build_valid_fvar_matrix<F>(template_vector));
          d.push_back(build_valid_fvar_matrix<F>(template_vector));
        }
      vector<VectorXvar> fd = F::template apply<vector<VectorXvar> >(d);

      EXPECT_EQ(d.size(), fd.size());
      EXPECT_EQ(d[i].size(), fd[i].size());
      test_fvar(F::apply_base(c[i](j)), c[i](j), fd[i](j), d[i](j));
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
    RowVectorXvar a = build_valid_fvar_matrix<F>(template_vector, i);
    RowVectorXvar b = build_valid_fvar_matrix<F>(template_vector, i);
    RowVectorXvar fb = F::template apply<RowVectorXvar>(b);
    EXPECT_EQ(b.size(), fb.size());
    test_fvar(F::apply_base(a(i)), a(i), fb(i), b(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<RowVectorXvar> c;
      vector<RowVectorXvar> d;
      for (size_t k = 0; k < vector_vector_size; ++k)
        if (k == i) { 
          c.push_back(build_valid_fvar_matrix<F>(template_vector, j));
          d.push_back(build_valid_fvar_matrix<F>(template_vector, j));
        }
        else {
          c.push_back(build_valid_fvar_matrix<F>(template_vector));
          d.push_back(build_valid_fvar_matrix<F>(template_vector));
        }
      vector<RowVectorXvar> fd = 
          F::template apply<vector<RowVectorXvar> >(d);

      EXPECT_EQ(d.size(), fd.size());
      EXPECT_EQ(d[i].size(), fd[i].size());
      test_fvar(F::apply_base(c[i](j)), c[i](j), fd[i](j), d[i](j));
    }
  }
}

template <typename F>
void expect_values() {
  using stan::math::fvar;
  using stan::math::var;
  using std::vector;

  expect_scalar_value<F, fvar<var> >();
  expect_scalar_value<F, fvar<fvar<var> > >();
  expect_std_vector_values<F, fvar<var> >();
  expect_std_vector_values<F, fvar<fvar<var> > >();
  expect_matrix_values<F, fvar<var> >();
  expect_matrix_values<F, fvar<fvar<var> > >();
  expect_vector_values<F, fvar<var> >();
  expect_vector_values<F, fvar<fvar<var> > >();
  expect_row_vector_values<F, fvar<var> >();
  expect_row_vector_values<F, fvar<fvar<var> > >();
}

#endif
