#ifndef TEST_UNIT_MATH_OPENCL_KERNEL_GENERATOR_REFERENCE_KERNEL_HPP
#define TEST_UNIT_MATH_OPENCL_KERNEL_GENERATOR_REFERENCE_KERNEL_HPP
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>

namespace stan {
namespace test {

/**
 * Loads file contents. File is assumed to be in
 * test/unit/math/opencl/kernel_generator/reference_kernels/
 * If file does not exist returns empty string.
 * @param filename Name of the file
 * @return file content
 */
inline std::string load_reference_kernel(const std::string& filename) {
  // check if we are in top level or test directory
  auto pwd = std::filesystem::current_path();
  while (pwd.stem() != "math" && pwd.has_parent_path()) {
    pwd = pwd.parent_path();
  }
  pwd /= "test/unit/math/opencl/kernel_generator/reference_kernels/" + filename;
  std::ifstream input(pwd);
  std::stringstream stream;
  stream << input.rdbuf();
  return stream.str();
}

/**
 * If \c STAN_TEST_KERNEL_GENERATOR_STORE_REFERENCE_KERNELS is defined saves the
 * kernel into a file with given name. Otherwise does nothing. File is assumed
 * to be in test/unit/math/opencl/kernel_generator/reference_kernels/
 * @param filename name of the file
 * @param kernel content to write in the file
 * @throw ios_base::failure File could not be opened or written.
 */
inline void store_reference_kernel_if_needed(const std::string& filename,
                                             const std::string& kernel) {
#ifdef STAN_TEST_KERNEL_GENERATOR_STORE_REFERENCE_KERNELS
  // check if we are in top level or test directory
  auto pwd = std::filesystem::current_path();
  while (pwd.stem() != "math" && pwd.has_parent_path()) {
    pwd = pwd.parent_path();
  }
  pwd /= "test/unit/math/opencl/kernel_generator/reference_kernels/" + filename;
  std::ofstream output;
  output.exceptions(~std::ofstream::goodbit);
  output.open(path);
  output << kernel;
#endif
}

}  // namespace test
}  // namespace stan

#endif  // TEST_UNIT_MATH_OPENCL_KERNEL_GENERATOR_REFERENCE_KERNEL_HPP
