#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>


namespace {
  /**
   * r: [1]
   * v: [2]
   * m: [7, 4]
   * ar: [2, 1]
   * aav: [2, 7, 7]
   * aam: [7, 4, 7, 4]
   * aar: [7, 4, 1]
   * aaar: [7, 4, 4, 1]
   */
  constexpr std::string_view test_json(
  "[{\"name\":\"real_p\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},"
  "{\"name\":\"vector_p\",\"type\":{\"name\":\"vector\",\"length\": 2 },\"block\":\"parameters\"},"
  "{\"name\":\"matrix_p\",\"type\":{\"name\":\"matrix\",\"rows\": 7 ,\"cols\": 4 },\"block\":\"parameters\"},"
  "{\"name\":\"arr_real_p\",\"type\":{\"name\":\"array\",\"length\": 2 ,\"element_type\":"
      "{\"name\":\"real\"}},\"block\":\"parameters\"},"
  "{\"name\":\"arr_vec_p\",\"type\":{\"name\":\"array\",\"length\": 2 ,\"element_type\":"
    "{\"name\":\"array\",\"length\": 7 ,\"element_type\":"
      "{\"name\":\"vector\",\"length\": 7 }}},\"block\":\"parameters\"},"
  "{\"name\":\"arr_arr_mat_p\",\"type\":{\"name\":\"array\",\"length\": 7 ,\"element_type\":"
    "{\"name\":\"array\",\"length\": 4 ,\"element_type\":"
      "{\"name\":\"matrix\",\"rows\": 7 ,\"cols\": 4 }}},\"block\":\"parameters\"},"
  "{\"name\":\"arr_arr_real_p\",\"type\":{\"name\":\"array\",\"length\": 7 ,\"element_type\":"
    "{\"name\":\"array\",\"length\": 4 ,\"element_type\":"
      "{\"name\":\"real\"}}},\"block\":\"parameters\"}]"
  "{\"name\":\"arr_arr_arr_real_p\",\"type\":{\"name\":\"array\",\"length\": 7 ,\"element_type\":"
    "{\"name\":\"array\",\"length\": 4 ,\"element_type\":"
      "{\"name\":\"array\",\"length\": 4 ,\"element_type\":"
       "{\"name\":\"real\"}}}},\"block\":\"parameters\"}]"
  "{\"name\":\"real_tp\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},"
  "{\"name\":\"vector_tp\",\"type\":{\"name\":\"vector\",\"length\": 2 },\"block\":\"transformed_parameters\"},"
  "{\"name\":\"matrix_tp\",\"type\":{\"name\":\"matrix\",\"rows\": 7 ,\"cols\": 4 },\"block\":\"transformed_parameters\"},"
  "{\"name\":\"arr_real_tp\",\"type\":{\"name\":\"array\",\"length\": 2 ,\"element_type\":"
  "{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},"
  "{\"name\":\"arr_vec_tp\",\"type\":{\"name\":\"array\",\"length\": 2 ,\"element_type\":"
    "{\"name\":\"array\",\"length\": 7 ,\"element_type\":"
      "{\"name\":\"vector\",\"length\": 7 }}},\"block\":\"transformed_parameters\"},"
  "{\"name\":\"arr_arr_mat_tp\",\"type\":{\"name\":\"array\",\"length\": 7 ,\"element_type\":"
    "{\"name\":\"array\",\"length\": 4 ,\"element_type\":"
      "{\"name\":\"matrix\",\"rows\": 7 ,\"cols\": 4 }}},\"block\":\"transformed_parameters\"},"
  "{\"name\":\"real_gq\",\"type\":{\"name\":\"real\"},\"block\":\"generated_quantities\"},"
  "{\"name\":\"vector_gq\",\"type\":{\"name\":\"vector\",\"length\": 2 },\"block\":\"generated_quantities\"},"
  "{\"name\":\"matrix_gq\",\"type\":{\"name\":\"matrix\",\"rows\": 7 ,\"cols\": 4 },\"block\":\"generated_quantities\"},"
  "{\"name\":\"arr_real_gq\",\"type\":{\"name\":\"array\",\"length\": 2 ,\"element_type\":"
    "{\"name\":\"real\"}},\"block\":\"generated_quantities\"},"
  "{\"name\":\"arr_vec_gq\",\"type\":{\"name\":\"array\",\"length\": 2 ,\"element_type\":"
    "{\"name\":\"array\",\"length\": 7 ,\"element_type\":"
      "{\"name\":\"vector\",\"length\": 7 }}},\"block\":\"generated_quantities\"},"
  "{\"name\":\"arr_arr_mat_gq\",\"type\":{\"name\":\"array\",\"length\": 7 ,\"element_type\":"
    "{\"name\":\"array\",\"length\": 4 ,\"element_type\":"
      "{\"name\":\"matrix\",\"rows\": 7 ,\"cols\": 4 }}},\"block\":\"generated_quantities\"}]");

TEST(OpenCLDeserializer, json_parse) {
  auto param_info = stan::math::internal::extract_parameter_types(test_json);
}

TEST(OpenCLDeserializer, deserializer) {
  const auto alignment_in_doubles = (stan::math::opencl_context.device()[0].getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() / 8)/sizeof(double);
  using stan::math::var_value;
  using stan::math::matrix_cl;
  using stan::math::var;
  using varmat_cl = var_value<matrix_cl<double>>;
  using std_vec_varmat_cl = std::vector<varmat_cl>;
  auto param_info = stan::math::internal::extract_parameter_types(test_json);
  std::size_t total_size = 0;
  for (auto& dim_i : param_info) {
    std::size_t dim_size = 1;
    for (auto dim_j : dim_i.second) {
      dim_size *= dim_j;
    }
    total_size += dim_size;
  }
  Eigen::VectorXd values = Eigen::VectorXd::LinSpaced(total_size, 1, total_size);
  varmat_cl deserialize_data = stan::math::make_deserializer_data(param_info, values);
  Eigen::MatrixXd from_data = stan::math::from_matrix_cl(deserialize_data.val());
  stan::math::deserializer_cl<varmat_cl> des_test(deserialize_data);
  Eigen::Index input_idx = 0;
  auto scalar_1 = des_test.template read<var>(1);
  EXPECT_EQ(scalar_1.val(), values(input_idx++));
  auto vec_1 = des_test.template read<varmat_cl>(2);
  auto host_vec_1 = stan::math::from_matrix_cl(vec_1.val());
  for (Eigen::Index i = 0; i < 2; i++) {
    EXPECT_EQ(host_vec_1(i), values(input_idx++));
  }
  auto mat_1 = des_test.template read<varmat_cl>(7, 4);
  auto host_mat_1 = stan::math::from_matrix_cl(mat_1.val());
  for (Eigen::Index j = 0; j < 4; j++) {
    for (Eigen::Index i = 0; i < 7; i++) {
      EXPECT_EQ(host_mat_1(i, j), values(input_idx++));
    }
  }
  auto std_vec_to_vec = des_test.template read<varmat_cl>(2);
  auto host_std_vec_to_vec = stan::math::from_matrix_cl(std_vec_to_vec.val());
  for (Eigen::Index i = 0; i < 2; i++) {
    EXPECT_EQ(host_std_vec_to_vec(i), values(input_idx++));
  }
  auto std_vec_std_vec_vec_1 = des_test.template read<std::vector<std_vec_varmat_cl>>(2, 7, 7);
  std::vector<std::vector<Eigen::VectorXd>> host_std_vec_std_vec_vec_1;
  for (Eigen::Index i = 0; i < 2; i++) {
    std::vector<Eigen::VectorXd> x;
    x.reserve(7);
    for (Eigen::Index j = 0; j < 7; j++) {
      x.push_back(stan::math::from_matrix_cl(std_vec_std_vec_vec_1[i][j].val()));
    }
    host_std_vec_std_vec_vec_1.push_back(std::move(x));

  }
  for (Eigen::Index k = 0; k < 2; k++) {
    for (Eigen::Index j = 0; j < 7; j++) {
      for (Eigen::Index i = 0; i < 7; i++) {
        EXPECT_EQ(host_std_vec_std_vec_vec_1[k][j](i), values(input_idx++));
      }
    }
  }
  auto std_vec_std_vec_mat_2 = des_test.template read<std::vector<std_vec_varmat_cl>>(7, 4, 7, 4);
  std::vector<std::vector<Eigen::MatrixXd>> host_std_vec_std_vec_mat_2;
  for (Eigen::Index i = 0; i < 7; i++) {
    std::vector<Eigen::MatrixXd> x;
    x.reserve(7);
    for (Eigen::Index j = 0; j < 4; j++) {
      x.push_back(stan::math::from_matrix_cl(std_vec_std_vec_mat_2[i][j].val()));
    }
    host_std_vec_std_vec_mat_2.push_back(std::move(x));
  }
  for (Eigen::Index k = 0; k < 7; k++) {
    for (Eigen::Index j = 0; j < 4; j++) {
      for (Eigen::Index l = 0; l < 4; l++) {
        for (Eigen::Index i = 0; i < 7; i++) {
          EXPECT_EQ(host_std_vec_std_vec_mat_2[k][j](i, l), values(input_idx++));
        }
      }
    }
  }
  auto std_vec_std_vec_to_mat = des_test.template read<varmat_cl>(7, 4);
  Eigen::MatrixXd host_std_vec_std_vec_to_mat = stan::math::from_matrix_cl(std_vec_std_vec_to_mat.val());
  for (Eigen::Index j = 0; j < 4; j++) {
    for (Eigen::Index i = 0; i < 7; i++) {
      EXPECT_EQ(host_std_vec_std_vec_to_mat(i, j), values(input_idx++));
    }
  }
  auto std_vec_std_vec_stdvec_to_mat = des_test.template read<std::vector<varmat_cl>>(7, 4, 4);
  std::vector<Eigen::MatrixXd> host_std_vec_std_vec_stdvec_to_mat;
  for (Eigen::Index i = 0; i < 7; i++) {
    host_std_vec_std_vec_stdvec_to_mat.push_back(stan::math::from_matrix_cl(std_vec_std_vec_stdvec_to_mat[i].val()));
  }
  for (Eigen::Index k = 0; k < 7; k++) {
    for (Eigen::Index j = 0; j < 4; j++) {
      for (Eigen::Index i = 0; i < 4; i++) {
        EXPECT_EQ(host_std_vec_std_vec_stdvec_to_mat[k](i, j), values(input_idx++));
      }
    }
  }
  deserialize_data.adj() += 2.0 * deserialize_data.val();
  Eigen::VectorXd unpadded_host_data = Eigen::VectorXd::Zero(values.size());
  extract_deserializer_data(unpadded_host_data, param_info, deserialize_data);
  for (Eigen::Index i = 0; i < values.size(); i++) {
    EXPECT_EQ(unpadded_host_data(i), 2.0 * values(i));
  }
}


template <typename T>
struct deserializer {
  std::reference_wrapper<T> vals_;
  Eigen::Index position_;
  explicit deserializer(T& v_vals)
      : position_(0), vals_(v_vals) {}
  template <typename U, stan::require_stan_scalar_t<U>* = nullptr>
  auto read() {
    auto& vals = vals_.get();
    auto res = vals(position_);
    position_++;
    return res;
  }
  template <typename U, stan::require_eigen_t<U>* = nullptr>
  auto read(Eigen::Index r, Eigen::Index c) {
    auto& vals = vals_.get();
    std::decay_t<U> res(r, c);
    for (Eigen::Index j = 0; j < c; j++) {
      for (Eigen::Index i = 0; i < r; i++) {
        res(i, j) = vals(position_);
        position_++;
      }
    }
    return res;
  }
  template <typename U, stan::require_eigen_t<U>* = nullptr>
  auto read(Eigen::Index r) {
    auto& vals = vals_.get();
    std::decay_t<U> res(r);
    for (Eigen::Index i = 0; i < r; i++) {
      res(i) = vals(position_);
      position_++;
    }
    return res;
  }
  template <typename U, typename... Idxs, stan::require_std_vector_t<U>* = nullptr>
  auto read(Eigen::Index r, Idxs... dims) {
    auto& vals = vals_.get();
    U res;
    res.reserve(r);
    for (Eigen::Index i = 0; i < r; i++) {
      res.push_back(this->template read<stan::value_type_t<U>>(dims...));
    }
    return res;
  }
};
TEST(OpenCLDeserializer, grad_test) {
  const auto alignment_in_doubles = (stan::math::opencl_context.device()[0].getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() / 8)/sizeof(double);
  using stan::math::var_value;
  using stan::math::matrix_cl;
  using stan::math::var;
  using varmat_cl = var_value<matrix_cl<double>>;
  using std_vec_varmat_cl = std::vector<varmat_cl>;
  auto param_info = stan::math::internal::extract_parameter_types(test_json);
  std::size_t total_size = 0;
  for (auto& dim_i : param_info) {
    std::size_t dim_size = 1;
    for (auto dim_j : dim_i.second) {
      dim_size *= dim_j;
    }
    total_size += dim_size;
  }
  Eigen::VectorXd values = Eigen::VectorXd::LinSpaced(total_size, 1, total_size);
  varmat_cl deserialize_data = stan::math::make_deserializer_data(param_info, values);
  Eigen::MatrixXd from_data = stan::math::from_matrix_cl(deserialize_data.val());
  using stan::math::sum;
  using stan::math::add;
  using stan::math::sin;
  using stan::math::multiply;
  auto grad_f = [&](auto&& x) {
    constexpr bool is_var_cl = std::is_same<std::decay_t<decltype(x)>, varmat_cl>::value;
    using matrix_t = std::conditional_t<is_var_cl, var_value<stan::math::matrix_cl<double>>, Eigen::Matrix<var, -1, -1>>;
    using vector_t = std::conditional_t<is_var_cl, var_value<stan::math::matrix_cl<double>>, Eigen::Matrix<var, -1, 1>>;
    using deserializer_t = std::conditional_t<is_var_cl, stan::math::deserializer_cl<varmat_cl>, deserializer<Eigen::Matrix<var, -1, 1>>>;
    using std_vec_matrix_t = std::vector<matrix_t>;
    using std_vec_vec_t = std::vector<vector_t>;
    deserializer_t des_test(x);
    auto scalar_1 = des_test.template read<var>();
    auto vec_1 = des_test.template read<vector_t>(2);
    auto mat_1 = des_test.template read<matrix_t>(7, 4);
    auto std_vec_to_vec = des_test.template read<vector_t>(2);
    auto std_vec_std_vec_vec_1 = des_test.template read<std::vector<std_vec_vec_t>>(2, 7, 7);
    auto std_vec_std_vec_mat_2 = des_test.template read<std::vector<std_vec_matrix_t>>(7, 4, 7, 4);
    auto std_vec_std_vec_to_mat = des_test.template read<matrix_t>(7, 4);
    auto std_vec_std_vec_stdvec_to_mat = des_test.template read<std_vec_matrix_t>(7, 4, 4);
    auto res = scalar_1 + sum(add(vec_1, std_vec_to_vec));
    for (int i = 0; i < 7; i++) {
      res += sum(
        add(
          add(
            add(multiply(std_vec_std_vec_mat_2[i][0], std_vec_std_vec_stdvec_to_mat[i]),
              multiply(std_vec_std_vec_mat_2[i][1], std_vec_std_vec_stdvec_to_mat[i])),
            multiply(std_vec_std_vec_mat_2[i][2], std_vec_std_vec_stdvec_to_mat[i])),
          multiply(std_vec_std_vec_mat_2[i][3], std_vec_std_vec_stdvec_to_mat[i])));
    }
    return res;
  };
  auto res = grad_f(deserialize_data);
  stan::math::grad(res.vi_);
  Eigen::VectorXd unpadded_host_adj = Eigen::VectorXd::Zero(values.size());
  extract_deserializer_data(unpadded_host_adj, param_info, deserialize_data);
  stan::math::recover_memory();
  Eigen::Matrix<var, -1, 1> host_matrix = values;
  auto host_res = grad_f(host_matrix);
  stan::math::grad(host_res.vi_);
  Eigen::VectorXd host_adj = host_matrix.adj();
  stan::math::recover_memory();
  for (Eigen::Index i = 0; i < host_adj.size(); i++) {
    EXPECT_EQ(host_adj(i), unpadded_host_adj(i));
  }
}

TEST(OpenCLDeserializer, gradient_fun_test) {
  const auto alignment_in_doubles = (stan::math::opencl_context.device()[0].getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() / 8)/sizeof(double);
  using stan::math::var_value;
  using stan::math::matrix_cl;
  using stan::math::var;
  using varmat_cl = var_value<matrix_cl<double>>;
  using std_vec_varmat_cl = std::vector<varmat_cl>;
  auto param_info = stan::math::internal::extract_parameter_types(test_json);
  std::size_t total_size = 0;
  for (auto& dim_i : param_info) {
    std::size_t dim_size = 1;
    for (auto dim_j : dim_i.second) {
      dim_size *= dim_j;
    }
    total_size += dim_size;
  }
  Eigen::VectorXd values = Eigen::VectorXd::LinSpaced(total_size, 1, total_size);
  matrix_cl<double> deserialize_data = stan::math::make_deserializer_data(param_info, values);
  using stan::math::sum;
  using stan::math::add;
  using stan::math::sin;
  using stan::math::multiply;
  auto grad_f = [&](auto&& x) {
    constexpr bool is_var_cl = std::is_same<std::decay_t<decltype(x)>, varmat_cl>::value;
    using matrix_t = std::conditional_t<is_var_cl, var_value<stan::math::matrix_cl<double>>, Eigen::Matrix<var, -1, -1>>;
    using vector_t = std::conditional_t<is_var_cl, var_value<stan::math::matrix_cl<double>>, Eigen::Matrix<var, -1, 1>>;
    using deserializer_t = std::conditional_t<is_var_cl, stan::math::deserializer_cl<varmat_cl>, deserializer<Eigen::Matrix<var, -1, 1>>>;
    using std_vec_matrix_t = std::vector<matrix_t>;
    using std_vec_vec_t = std::vector<vector_t>;
    deserializer_t des_test(x);
    auto scalar_1 = des_test.template read<var>();
    auto vec_1 = des_test.template read<vector_t>(2);
    auto mat_1 = des_test.template read<matrix_t>(7, 4);
    auto std_vec_to_vec = des_test.template read<vector_t>(2);
    auto std_vec_std_vec_vec_1 = des_test.template read<std::vector<std_vec_vec_t>>(2, 7, 7);
    auto std_vec_std_vec_mat_2 = des_test.template read<std::vector<std_vec_matrix_t>>(7, 4, 7, 4);
    auto std_vec_std_vec_to_mat = des_test.template read<matrix_t>(7, 4);
    auto std_vec_std_vec_stdvec_to_mat = des_test.template read<std_vec_matrix_t>(7, 4, 4);
    auto res = scalar_1 + sum(add(vec_1, std_vec_to_vec));
    for (int i = 0; i < 7; i++) {
      res += sum(
        add(
          add(
            add(multiply(std_vec_std_vec_mat_2[i][0], std_vec_std_vec_stdvec_to_mat[i]),
              multiply(std_vec_std_vec_mat_2[i][1], std_vec_std_vec_stdvec_to_mat[i])),
            multiply(std_vec_std_vec_mat_2[i][2], std_vec_std_vec_stdvec_to_mat[i])),
          multiply(std_vec_std_vec_mat_2[i][3], std_vec_std_vec_stdvec_to_mat[i])));
    }
    return res;
  };
  double res_val{0};
  matrix_cl<double> grad_val(deserialize_data.rows(), deserialize_data.cols());
  stan::math::gradient(grad_f, deserialize_data, res_val, grad_val);
  Eigen::VectorXd unpadded_host_adj = Eigen::VectorXd::Zero(values.size());
  extract_deserializer_data(unpadded_host_adj, param_info, grad_val);
  stan::math::recover_memory();
  Eigen::Matrix<var, -1, 1> host_matrix = values;
  Eigen::VectorXd host_adj = Eigen::VectorXd::Zero(values.size());
  res_val = 0;
  stan::math::gradient(grad_f, values, res_val, host_adj);
  stan::math::recover_memory();
  for (Eigen::Index i = 0; i < host_adj.size(); i++) {
    EXPECT_EQ(host_adj(i), unpadded_host_adj(i));
  }
}
}
#endif
