#ifndef STAN_MATH_OPENCL_REV_DESERIALIZER_HPP
#define STAN_MATH_OPENCL_REV_DESERIALIZER_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan::math {
namespace internal {
// Helper: trim whitespace from both ends.
template <typename T>
inline auto trim(T&& s) {
  const char* whitespace = " \t\n\r";
  size_t start = s.find_first_not_of(whitespace);
  if (start == std::string::npos)
    return decltype(s)("");
  size_t end = s.find_last_not_of(whitespace);
  return s.substr(start, end - start + 1);
}

// Helper: Given a JSON string and a field name, extract the field’s string value.
template <typename T1, typename T2>
inline auto get_json_string_field(T1&& json, T2&& field) {
  std::string key = "\"" + std::string(field) + "\"";
  size_t pos = json.find(key);
  if (pos == std::string::npos)
    throw std::runtime_error("Field not found: " + std::string(field));
  pos = json.find(":", pos);
  if (pos == std::string::npos)
    throw std::runtime_error("Colon not found for field: " + std::string(field));
  pos = json.find("\"", pos);
  if (pos == std::string::npos)
    throw std::runtime_error("Opening quote not found for field: " + std::string(field));
  size_t start = pos + 1;
  size_t end = json.find("\"", start);
  if (end == std::string::npos)
    throw std::runtime_error("Closing quote not found for field: " + std::string(field));
  return json.substr(start, end - start);
}

// Helper: Given a JSON string and a field name, extract the field’s integer value.
template <typename T1, typename T2>
inline int get_json_int_field(T1&& json, T2&& field) {
  std::string key = "\"" + std::string(field) + "\"";
  size_t pos = json.find(key);
  if (pos == std::string::npos)
    throw std::runtime_error("Field not found: " + std::string(field));
  pos = json.find(":", pos);
  if (pos == std::string::npos)
    throw std::runtime_error("Colon not found for field: " + std::string(field));
  pos++; // move past the colon
  // Skip whitespace
  while (pos < json.size() && std::isspace(json[pos]))
    pos++;
  size_t start = pos;
  while (pos < json.size() && (std::isdigit(json[pos]) || json[pos]=='-' || json[pos]=='+'))
    pos++;
  if (start == pos)
    throw std::runtime_error("No numeric value for field: " + std::string(field));
  auto num_str = json.substr(start, pos - start);
  return std::stoi(num_str);
}

// Helper: Given a string and a starting position (pointing at a '{'),
// extract the entire JSON object (using brace matching).
template <typename T>
inline auto extract_json_object(T&& s, size_t pos) {
  if (s[pos] != '{')
    throw std::runtime_error("Expected '{' at position " + std::to_string(pos));
  int count = 0;
  size_t start = pos;
  size_t i = pos;
  for (; i < s.size(); ++i) {
    if (s[i] == '{') count++;
    else if (s[i] == '}') {
      count--;
      if (count == 0)
        break;
    }
  }
  if (count != 0)
    throw std::runtime_error("Unmatched braces in JSON object.");
  return s.substr(start, i - start + 1);
}

// Recursive function to parse a type JSON string and return a pair:
//   first: a short string key (r, v, m, a...)
//   second: a vector of ints representing the dimensions.
template <typename T>
std::pair<std::string, std::vector<int>> parse_type_string(T&& s) {
  std::pair<std::string, std::vector<int>> result;
  auto name = get_json_string_field(s, "name");
  if (name == "real") {
    result.first = "r";
    result.second = {1};
  } else if (name == "vector") {
    result.first = "v";
    int len = get_json_int_field(s, "length");
    result.second = {len};
  } else if (name == "matrix") {
    result.first = "m";
    int rows = get_json_int_field(s, "rows");
    int cols = get_json_int_field(s, "cols");
    result.second = {rows, cols};
  } else if (name == "array") {
    // For an array, get its "length" and then recursively parse its element_type.
    int len = get_json_int_field(s, "length");
    size_t pos = s.find("\"element_type\"");
    if (pos == std::string::npos)
      throw std::runtime_error("array missing element_type");
    pos = s.find(":", pos);
    if (pos == std::string::npos)
      throw std::runtime_error("array missing colon after element_type");
    pos = s.find("{", pos);
    if (pos == std::string::npos)
      throw std::runtime_error("array missing '{' for element_type");
    auto elem_type_str = extract_json_object(s, pos);
    auto inner = parse_type_string(elem_type_str);
    // Prepend an "a" to indicate array.
    result.first = "a" + inner.first;
    // The dimensions: first the array length, then the inner dimensions.
    result.second.push_back(len);
    for (auto val : inner.second)
      result.second.push_back(val);
  } else {
    throw std::runtime_error("Unknown type name: " + name);
  }
  return result;
}

// A minimal parser that scans the JSON string for top-level objects (using brace matching)
// and extracts the "type" field from those objects if their "block" field is "parameters"
// or "transformed_parameters". It returns a vector of the raw JSON for the type field.
template <typename T>
std::vector<std::string> parse_parameter_types(T&& json) {
  std::vector<std::string> raw_types;
  size_t pos = 0;
  while ((pos = json.find("{", pos)) != std::string::npos) {
    size_t start = pos;
    int braceCount = 0;
    size_t end = pos;
    for (; end < json.size(); ++end) {
      if (json[end] == '{') ++braceCount;
      else if (json[end] == '}') {
        --braceCount;
        if (braceCount == 0)
          break;
      }
    }
    if (braceCount != 0)
      throw std::runtime_error("Malformed JSON: unmatched braces.");
    auto objectStr = json.substr(start, end - start + 1);
    pos = end + 1;
    // Only process objects whose "block" is "parameters" or "transformed_parameters"
    if (objectStr.find("\"block\":\"parameters\"") == std::string::npos)
      continue;
    // Find the "type" field.
    size_t typePos = objectStr.find("\"type\":");
    if (typePos == std::string::npos)
      continue;
    size_t typeObjStart = objectStr.find("{", typePos);
    if (typeObjStart == std::string::npos)
      continue;
    auto typeStr = extract_json_object(objectStr, typeObjStart);
    raw_types.push_back(std::string(trim(typeStr)));
  }
  return raw_types;
}

/**
 * Given a JSON string (as produced by get_constrained_sizedtypes), parse it and return
 * a vector of pairs. Each pair consists of a short type code and a vector of integers describing
 * the dimensions. For example, for the sample JSON the return value is:
 *
 *   { {"r", {1}},
 *     {"v", {2}},
 *     {"m", {7,4}},
 *     {"ar", {2}},
 *     {"aav", {2,7,7}},
 *     {"aam", {7,4,7,4}},
 *     {"r", {1}},
 *     {"v", {2}},
 *     {"m", {7,4}},
 *     {"ar", {2}},
 *     {"aav", {2,7,7}},
 *     {"aam", {7,4,7,4}} }
 *
 * @param json The JSON string to parse.
 * @return A vector of pairs of type (short type code, dimensions vector).
 */
template <typename T>
std::vector<std::pair<std::string, std::vector<int>>> extract_parameter_types(T&& json) {
  std::vector<std::pair<std::string, std::vector<int>>> ret;
  auto raw_types = parse_parameter_types(json);
  for (const auto & s : raw_types) {
    ret.push_back(parse_type_string(s));
  }
  return ret;
}
}

// Example function to compute the padded offset.
// alignment_in_doubles is the alignment requirement in units of doubles.
inline size_t compute_padded_offset(size_t current_offset, size_t alignment_in_doubles) {
  size_t remainder = current_offset % alignment_in_doubles;
  if (remainder == 0)
    return current_offset;
  else
    return current_offset + (alignment_in_doubles - remainder);
}

// TODO(Steve): iterate over subsets better here.
/**
r: [1, ]
v: [2, ]
m: [7, 4, ]
ar: [2, 1, ]
aav: [2, 7, 7, ]
aam: [7, 4, 7, 4, ]

 */
/**
 * Create a `var_value<matrix_cl<double>>` with the correct padding and values
 * for the subbuffers.
 *
 * @tparam Model The model class which has methods `get_param_names` and `get_dims`.
 * @param model A pointer to a Stan model.
 * @param values The Eigen vector whose values will be placed in the correct positions
 * for the values of the `var_value<matrix_cl<double>>`
 * @return a var_value<matrix_cl<double>> whose buffer is allocated with proper padding.
 */
template <typename Dims>
inline auto make_deserializer_data(const Dims& dims, const Eigen::VectorXd& values) {

  /*
   * `CL_DEVICE_MEM_BASE_ADDR_ALIGN` is the alignment requirement (in bits) for
   * sub-buffer offsets. The minimum value
   * is the size (in bits) of the largest OpenCL built-in data type supported
   * by the device (long16 in FULL profile, long16 or int16 in EMBEDDED profile)
   * for devices that are not of type CL_DEVICE_TYPE_CUSTOM.
   * The standard range for this is from 128 to 4096 bits (16 to 512 bytes).
   */
  const size_t alignment_in_doubles = (opencl_context.device()[0].getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() / 8)/sizeof(double);

  // STEP 3. Compute padded offsets for each parameter.
  size_t total_doubles = 0;
  std::vector<size_t> offsets;  // will hold the starting offset (in doubles) for each parameter
  for (const auto& dim : dims) {
    // Compute the number of elements for this parameter.
    Eigen::Index param_size = 1;
    std::cout << "dim.first: " << dim.first << ":\n";
    Eigen::Index i = dim.second.size() - 1;
    for (; i >= 0; --i) {
      if (dim.first[i] == 'a') {
        break;
      }
      param_size *= dim.second[i];
    }
    // Will only happen if the dims contain 'a'.
    std::cout << "LINE: " << __LINE__ << std::endl;
    Eigen::Index num_arrays = 1;
    for (; i >= 0; --i) {
      num_arrays *= dim.second[i];
    }
    // for array[N] real we treat it as one array of size N. Otherwise we need to pad each one which can get expensive.
    if (num_arrays > 1 && dim.first.back() == 'r') {
      param_size *= num_arrays;
      num_arrays = 1;
    }
    std::cout << "LINE: " << __LINE__ << std::endl;
    for (; num_arrays > 0; num_arrays--) {
      // Ensure the starting offset is aligned.
      total_doubles = compute_padded_offset(total_doubles, alignment_in_doubles);
      offsets.push_back(total_doubles);
      total_doubles += param_size;
    }
  }
  std::cout << "Offsets: [";
  for (auto offset : offsets) {
    std::cout << offset << ", ";
  }
  std::cout << "]\n";
  std::cout << "total_doubles: " << total_doubles << std::endl;
  // Now total_doubles is the size of our padded global buffer.

  // STEP 4. Allocate a padded host vector and copy in the values.
  Eigen::VectorXd padded_data(total_doubles);
  padded_data.setZero();  // initialize to zero so padding regions are defined
  Eigen::Index input_pos = 0;
  std::cout << "===========\n";
  for (size_t i = 0; i < dims.size(); ++i) {
    std::cout << "iter: " << i << std::endl;
    std::cout << "input_pos: " << input_pos << std::endl;
    size_t param_size = 1;
    std::cout << "dim(" << i << "): " << dims[i].first << ": [";
    for (size_t d : dims[i].second) {
      std::cout << d << ", ";
      param_size *= d;
    }
    std::cout << "]" << std::endl;
    std::cout << "param_size: " << param_size << std::endl;
    // Ensure we have enough input values.
    std::cout << "input_pos + param_size: " << (input_pos + param_size) << std::endl;
    assert(input_pos + param_size <= static_cast<size_t>(values.size()));
    // Copy the parameter's data from the input vector into the correct subregion.
    padded_data.segment(offsets[i], param_size) = values.segment(input_pos, param_size);
    input_pos += param_size;
    std::cout << "===============\n";
  }
  // Optionally, you might also want to check that input_pos == values.size()

  // STEP 5. Create a matrix_cl from the padded_data.
  // We assume that the matrix_cl constructor can accept an Eigen vector
  // and will allocate a buffer on the GPU (with SVM or via OpenCL buffers).
  stan::math::matrix_cl<double> gpu_buffer(padded_data);

  // STEP 6. Wrap the matrix_cl in a var_value and return it.
  return stan::math::var_value<stan::math::matrix_cl<double>>(std::move(gpu_buffer));
}

template <typename Mat>
struct deserializer_cl {
  std::reference_wrapper<Mat> mat_;
  std::size_t alignment_in_doubles_{64};
  std::size_t current_offset_{0};
  deserializer_cl(Mat& mat, std::size_t alignment_in_doubles)
      : mat_(mat), alignment_in_doubles_(alignment_in_doubles), current_offset_(0) {}

  explicit deserializer_cl(Mat& mat)
      : mat_(mat),
      alignment_in_doubles_((opencl_context.device()[0].getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() / 8)/sizeof(double)),
      current_offset_(0) {}
  inline void align_offset() {
    const size_t remainder = current_offset_ % alignment_in_doubles_;
    if (remainder == 0)
      return;
    else
      current_offset_ = current_offset_ + (alignment_in_doubles_ - remainder);
  }

  template <typename T, require_var_t<T>* = nullptr, require_matrix_cl_t<value_type_t<T>>* = nullptr>
  inline auto read(Eigen::Index amt) {
    align_offset();
    auto& mat = mat_.get();
    auto& val_buffer = mat.vi_->val_.buffer();
    cl_buffer_region region;
    region.origin = current_offset_ * sizeof(double);
    region.size = amt * sizeof(double);
    cl_int err = CL_SUCCESS;
    cl::Buffer val_sub_buf = val_buffer.createSubBuffer(
      CL_MEM_READ_WRITE,
      CL_BUFFER_CREATE_TYPE_REGION,
      &region,
      &err
    );
    if (err != CL_SUCCESS) {
      throw std::domain_error("Failed to create sub-buffer");
    }
    auto& adj_buffer = mat.vi_->adj_.buffer();
    cl::Buffer adj_sub_buf = adj_buffer.createSubBuffer(
      CL_MEM_READ_WRITE,
      CL_BUFFER_CREATE_TYPE_REGION,
      &region,
      &err
    );
    if (err != CL_SUCCESS) {
      throw std::domain_error("Failed to create sub-buffer");
    }
    matrix_cl<double> val(val_sub_buf, amt, 1);
    matrix_cl<double> adj(adj_sub_buf, amt, 1);
    current_offset_ += amt;
    align_offset();
    return var_value<matrix_cl<double>>(new vari_value<matrix_cl<double>>(std::move(val), std::move(adj)));
  }
  // TODO(Steve): This can be one variadic parameter pack
  template <typename T, require_var_t<T>* = nullptr, require_matrix_cl_t<value_type_t<T>>* = nullptr>
  inline auto read(Eigen::Index rows, Eigen::Index cols) {
    const auto amt = rows * cols;
    align_offset();
    auto& mat = mat_.get();
    auto& val_buffer = mat.vi_->val_.buffer();
    cl_buffer_region region;
    region.origin = current_offset_ * sizeof(double);
    region.size = amt * sizeof(double);
    cl_int err = CL_SUCCESS;
    cl::Buffer val_sub_buf = val_buffer.createSubBuffer(
      CL_MEM_READ_WRITE,
      CL_BUFFER_CREATE_TYPE_REGION,
      &region,
      &err
    );
    if (err != CL_SUCCESS) {
      throw std::domain_error("Failed to create sub-buffer");
    }
    auto& adj_buffer = mat.vi_->adj_.buffer();
    cl::Buffer adj_sub_buf = adj_buffer.createSubBuffer(
      CL_MEM_READ_WRITE,
      CL_BUFFER_CREATE_TYPE_REGION,
      &region,
      &err
    );
    if (err != CL_SUCCESS) {
      throw std::domain_error("Failed to create sub-buffer");
    }
    matrix_cl<double> val(std::move(val_sub_buf), rows, cols);
    matrix_cl<double> adj(std::move(adj_sub_buf), rows, cols);
    current_offset_ += amt;
    align_offset();
    return var_value<matrix_cl<double>>(new vari_value<matrix_cl<double>>(std::move(val), std::move(adj)));
  }

template <typename T, require_var_t<T>* = nullptr, require_arithmetic_t<value_type_t<T>>* = nullptr>
inline auto read(Eigen::Index amt) {
  if (amt != 1) {
    throw std::domain_error("INTERNAL ERROR: Scalar read must have amount of 1");
  }
  align_offset();
  // 1. Create a sub-buffer for a single double from the value buffer.
  auto& val_buffer = mat_.get().vi_->val_.buffer();
  cl_buffer_region region;
  region.origin = current_offset_ * sizeof(double);
  region.size   = sizeof(double);
  cl_int err = CL_SUCCESS;
  cl::Buffer val_sub_buf = val_buffer.createSubBuffer(
    CL_MEM_READ_WRITE,
    CL_BUFFER_CREATE_TYPE_REGION,
    &region,
    &err
  );
  if (err != CL_SUCCESS)
    throw std::domain_error("Failed to create sub-buffer for scalar value");

  // Read the single double from the device.
  double scalar_val = 0.0;
  stan::math::opencl_context.queue().enqueueReadBuffer(val_sub_buf, CL_TRUE, 0, sizeof(double), &scalar_val);

  // 2. Create a new vari object on the CPU using the scalar value.
  // Here we use Stan's standard vari (which holds a double value).
  // (If you need a GPU-friendly type, you may need to define your own vari_cl_scalar.)
  auto* var_ptr = new vari(scalar_val);
  var ret(var_ptr);

  // Capture the current offset before incrementing.
  size_t scalar_offset = current_offset_;

  // 3. Schedule a reverse pass callback that writes the CPU-computed adjoint
  // back to the corresponding position in the GPU's adjoint buffer.
  reverse_pass_callback([serial_mat = this->mat_,
                         scalar_offset,
                         var_ptr]() mutable {
    double grad = var_ptr->adj();
    auto& adj_buffer = serial_mat.get().vi_->adj_.buffer();
    cl_buffer_region region;
    region.origin = scalar_offset * sizeof(double);
    region.size   = sizeof(double);
    cl_int err = CL_SUCCESS;
    cl::Buffer adj_sub_buf = adj_buffer.createSubBuffer(
      CL_MEM_READ_WRITE,
      CL_BUFFER_CREATE_TYPE_REGION,
      &region,
      &err
    );
    if (err != CL_SUCCESS)
      throw std::domain_error("Failed to create sub-buffer for scalar adjoint");

    // Write the computed gradient back into the GPU adjoint memory.
    stan::math::opencl_context.queue().enqueueWriteBuffer(adj_sub_buf, CL_TRUE, 0, sizeof(double), &grad);
  });

  // 4. Advance the current offset and re-align.
  current_offset_ += 1;
  align_offset();
  return ret;
}
  template <typename T, typename... Idxs, require_std_vector_t<T>* = nullptr>
  inline auto read(Eigen::Index idx_1, Idxs... idxs) {
  align_offset();
    std::decay_t<T> vec;
    vec.reserve(idx_1);
    for (Eigen::Index i = 0; i < idx_1; ++i) {
      vec.emplace_back(this->read<value_type_t<std::decay_t<T>>>(idxs...));
    }
  align_offset();
    return vec;
  }

};

}
#endif
#endif
