#ifndef STAN_MATH_VULKAN_FROM_MATRIX_VK_HPP
#define STAN_MATH_VULKAN_FROM_MATRIX_VK_HPP
#ifdef STAN_VULKAN

#include <stan/math/vulkan/vulkan_context.hpp>
#include <stan/math/vulkan/matrix_vk.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> from_matrix_vk(const matrix_vk<T>& m) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(m.rows(), m.cols());
  T* data = static_cast<T*>(vulkan_context.device().mapMemory(m.buffer_memory(), 0, sizeof(T) * m.size()));
  std::copy(data, data + m.size(), result.data());
  return result;
}

}  // namespace math
}  // namespace stan

#endif
#endif
