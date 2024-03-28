#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/vulkan/vulkan_context.hpp>
#include <stan/math/vulkan/matrix_vk.hpp>
#include <stan/math/vulkan/vulkan_operation.hpp>
#include <stan/math/vulkan/from_matrix_vk.hpp>


int main() {
  stan::math::vulkan_context.printProperties();

  Eigen::MatrixXi data(3, 3);
  data << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  stan::math::matrix_vk<int> inp1(data);
  stan::math::matrix_vk<int> inp2(data);
  stan::math::matrix_vk<int> out(3, 3);

  stan::math::vulkan_operation add_op("stan/math/vulkan/shaders/add.spv");
  add_op.add_operands(out, inp1, inp2);
  add_op.eval_operation();

  Eigen::MatrixXi result = stan::math::from_matrix_vk(out);

  std::cout
    << "\nInput: \n" << data
    << "\nOutput: \n" << result
    << std::endl;

  return 0;
}
