#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
/*
TEST(OpenCLDeserializer, make_data) {
  std::vector<std::vector<size_t>> dims = {{2, 4}, {2, 2, 4}, {32}, {63}, {8, 4, 2}, {8, 4, 2, 2}};
  std::size_t total_size = 0;
  for (auto& dim_i : dims) {
    std::size_t dim_size = 1;
    for (auto dim_j : dim_i) {
      dim_size *= dim_j;
    }
    total_size += dim_size;
  }
  Eigen::VectorXd values = Eigen::VectorXd::LinSpaced(total_size, 1, 10);
  const auto alignment_in_doubles = (stan::math::opencl_context.device()[0].getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() / 8)/sizeof(double);
  size_t total_doubles = 0;
    std::vector<size_t> param_sizes;
  std::vector<size_t> offsets;  // will hold the starting offset (in doubles) for each parameter
  for (const auto& dim : dims) {
    // Compute the number of elements for this parameter.
    size_t param_size = 1;
    for (size_t d : dim) {
      param_size *= d;
    }
    // Ensure the starting offset is aligned.
    total_doubles = stan::math::compute_padded_offset(total_doubles, alignment_in_doubles);
    offsets.push_back(total_doubles);
    param_sizes.push_back(param_size);
    total_doubles += param_size;
  }
  std::cout << "Offsets: ";
  for (auto offset : offsets) {
    std::cout << offset << " ";
  }
  std::cout << std::endl;
  std::cout << "Alignment in doubles: " << alignment_in_doubles << std::endl;
  auto data = stan::math::make_deserializer_data(dims, values);
  EXPECT_EQ(data.val().rows(), total_doubles);
  EXPECT_EQ(data.val().cols(), 1);
  Eigen::VectorXd host_data = stan::math::from_matrix_cl(data.val());
  std::cout << "Host: " << host_data.size() << "\n" << host_data.transpose().eval() << std::endl;
  // Get the host buffer and attempt to make subbuffers from the offsets
    // Get the underlying OpenCL buffer from our matrix_cl.
  cl::Buffer& global_buffer = data.vi_->val_.buffer();
  cl::CommandQueue queue = stan::math::opencl_context.queue();

  // For each parameter, create a sub-buffer, read it back, and print its values.
  for (size_t i = 0; i < dims.size(); ++i) {
    size_t param_size = param_sizes[i];
    // Set up the region (offset and size in bytes)
    cl_buffer_region region;
    region.origin = offsets[i] * sizeof(double);
    region.size   = param_size * sizeof(double);

    cl_int err = CL_SUCCESS;
    cl::Buffer sub_buf = global_buffer.createSubBuffer(
      CL_MEM_READ_WRITE,
      CL_BUFFER_CREATE_TYPE_REGION,
      &region,
      &err
    );
    EXPECT_EQ(err, CL_SUCCESS);

    // Read the sub-buffer back to host.
    std::vector<double> sub_host(param_size, 0.0);
    err = queue.enqueueReadBuffer(sub_buf, CL_TRUE, 0, region.size, sub_host.data());
    EXPECT_EQ(err, CL_SUCCESS);

    // For verification, compare with the corresponding segment in the global padded host data.
    Eigen::VectorXd expected = host_data.segment(offsets[i], param_size);
    Eigen::Map<Eigen::VectorXd> sub_vec(sub_host.data(), sub_host.size());
    EXPECT_TRUE(sub_vec.isApprox(expected))
      << "Mismatch in parameter " << i;

    // Print the sub-buffer for visual inspection.
    std::cout << "Parameter " << i << " subbuffer (" << param_size << " values): "
              << sub_vec.transpose().eval() << std::endl;
  }
}
*/
namespace {
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
  std::cout << "param_info: \n" << std::endl;
  for (auto& x : param_info) {
    std::cout << x.first << ": [";
    for (auto& y : x.second) {
      std::cout << y << ", ";
    }
    std::cout << "]" << std::endl;
  }
  std::cout << std::endl;
}

TEST(OpenCLDeserializer, deserializer) {
  const auto alignment_in_doubles = (stan::math::opencl_context.device()[0].getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() / 8)/sizeof(double);
  std::cout << "Alignment in doubles: " << alignment_in_doubles << std::endl;
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
  Eigen::VectorXd values = Eigen::VectorXd::LinSpaced(total_size, 1, 10);
  for (auto& x : param_info) {
    std::cout << x.first << ": [";
    for (auto& y : x.second) {
      std::cout << y << ", ";
    }
    std::cout << "]" << std::endl;
  }
  std::cout << std::endl;
  varmat_cl deserialize_data = stan::math::make_deserializer_data(param_info, values);
  stan::math::deserializer_cl<varmat_cl> des_test(deserialize_data);
  std::cout << "Line: " << __LINE__ << std::endl;
  auto mat_11 = des_test.template read<var>(1);
  std::cout << "Line: " << __LINE__ << std::endl;
  auto std_vec_mat_1 = des_test.template read<varmat_cl>(2);
  std::cout << "Line: " << __LINE__ << std::endl;
  auto vec_1 = des_test.template read<varmat_cl>(7, 4);
  std::cout << "Line: " << __LINE__ << std::endl;
  auto vec_2 = des_test.template read<std::vector<var>>(2, 1);
  std::cout << "Line: " << __LINE__ << std::endl;
  auto std_vec_std_vec_vec_1 = des_test.template read<std::vector<std_vec_varmat_cl>>(2, 7, 7);
  std::cout << "Line: " << __LINE__ << std::endl;
  auto std_vec_std_vec_mat_2 = des_test.template read<std::vector<std_vec_varmat_cl>>(7, 4, 7, 4);
  std::cout << "Line: " << __LINE__ << std::endl;
/*
*/
}
}
#endif
