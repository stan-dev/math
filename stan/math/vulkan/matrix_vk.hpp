#ifndef STAN_MATH_VULKAN_MATRIX_VK_HPP
#define STAN_MATH_VULKAN_MATRIX_VK_HPP
#ifdef STAN_VULKAN

#include <stan/math/vulkan/vulkan_context.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename>
class matrix_vk;

template <typename T>
class matrix_vk {
 private:
  vk::Buffer buffer_vk_;
  vk::DeviceMemory buffer_memory_vk_;
  int rows_{0};
  int cols_{0};

 public:
  int rows() const { return rows_; }
  int cols() const { return cols_; }
  int size() const { return rows_ * cols_; }

  const vk::Buffer& buffer() const { return buffer_vk_; }
  const vk::DeviceMemory& buffer_memory() const { return buffer_memory_vk_; }

  matrix_vk() {};

  matrix_vk(int rows, int cols) : rows_(rows), cols_(cols) {
    vk::BufferCreateInfo buffer_info{};
    buffer_info.size = sizeof(T) * rows * cols;
    buffer_info.usage = vk::BufferUsageFlagBits::eStorageBuffer;
    buffer_info.sharingMode = vk::SharingMode::eExclusive;

    buffer_vk_ = vulkan_context.device().createBuffer(buffer_info);

    vk::MemoryRequirements mem_reqs = vulkan_context.device().getBufferMemoryRequirements(buffer_vk_);

    vk::MemoryAllocateInfo alloc_info{};
    alloc_info.allocationSize = mem_reqs.size;
    alloc_info.memoryTypeIndex = vulkan_context.memory_type_index();

    buffer_memory_vk_ = vulkan_context.device().allocateMemory(alloc_info);

    vulkan_context.device().bindBufferMemory(buffer_vk_, buffer_memory_vk_, 0);
  }

  matrix_vk(const matrix_vk<T>& m) : rows_(m.rows()), cols_(m.cols()) {
    vk::BufferCreateInfo buffer_info{};
    buffer_info.size = sizeof(T) * rows_ * cols_;
    buffer_info.usage = vk::BufferUsageFlagBits::eStorageBuffer;
    buffer_info.sharingMode = vk::SharingMode::eExclusive;

    buffer_vk_ = vulkan_context.device().createBuffer(buffer_info);

    vk::MemoryRequirements mem_reqs = vulkan_context.device().getBufferMemoryRequirements(buffer_vk_);

    vk::MemoryAllocateInfo alloc_info{};
    alloc_info.allocationSize = mem_reqs.size;
    alloc_info.memoryTypeIndex = vulkan_context.memory_type_index();

    buffer_memory_vk_ = vulkan_context.device().allocateMemory(alloc_info);

    vulkan_context.device().bindBufferMemory(buffer_vk_, buffer_memory_vk_, 0);

    vk::CommandBufferAllocateInfo alloc_info_cmd{};
    alloc_info_cmd.commandPool = vulkan_context.command_pool();
    alloc_info_cmd.level = vk::CommandBufferLevel::ePrimary;
    alloc_info_cmd.commandBufferCount = 1;

    vk::CommandBuffer cmd_buffer = vulkan_context.device().allocateCommandBuffers(alloc_info_cmd)[0];

    vk::CommandBufferBeginInfo begin_info{};
    begin_info.flags = vk::CommandBufferUsageFlagBits::eOneTimeSubmit;

    cmd_buffer.begin(begin_info);

    vk::BufferCopy copy_region{};
    copy_region.size = sizeof(T) * rows_ * cols_;

    cmd_buffer.copyBuffer(m.buffer(), buffer_vk_, copy_region);

    cmd_buffer.end();

    vk::SubmitInfo submit_info{};
    submit_info.commandBufferCount = 1;
    submit_info.pCommandBuffers = &cmd_buffer;

    vulkan_context.queue().submit(submit_info, nullptr);
    vulkan_context.queue().waitIdle();
  }

  matrix_vk(matrix_vk<T>&& m) : buffer_vk_(m.buffer_vk_), buffer_memory_vk_(m.buffer_memory_vk_), rows_(m.rows()), cols_(m.cols()) {
    m.buffer_vk_ = nullptr;
    m.buffer_memory_vk_ = nullptr;
    m.rows_ = 0;
    m.cols_ = 0;
  }

  matrix_vk(const std::vector<T>& vec) : rows_(vec.size()), cols_(1) {
    vk::BufferCreateInfo buffer_info{};
    buffer_info.size = sizeof(T) * rows_;
    buffer_info.usage = vk::BufferUsageFlagBits::eStorageBuffer;
    buffer_info.sharingMode = vk::SharingMode::eExclusive;

    buffer_vk_ = vulkan_context.device().createBuffer(buffer_info);

    vk::MemoryRequirements mem_reqs = vulkan_context.device().getBufferMemoryRequirements(buffer_vk_);

    vk::MemoryAllocateInfo alloc_info{};
    alloc_info.allocationSize = mem_reqs.size;
    alloc_info.memoryTypeIndex = vulkan_context.memory_type_index();

    buffer_memory_vk_ = vulkan_context.device().allocateMemory(alloc_info);

    vulkan_context.device().bindBufferMemory(buffer_vk_, buffer_memory_vk_, 0);

    T* data = static_cast<T*>(vulkan_context.device().mapMemory(buffer_memory_vk_, 0, sizeof(T) * rows_, vk::MemoryMapFlags()));
    std::copy(vec.begin(), vec.end(), data);
    vulkan_context.device().unmapMemory(buffer_memory_vk_);
  }

  template <int R, int C>
  matrix_vk(const Eigen::Matrix<T, R, C>& mat): rows_(mat.rows()), cols_(mat.cols()) {
    vk::BufferCreateInfo buffer_info{};
    buffer_info.size = sizeof(T) * rows_ * cols_;
    buffer_info.usage = vk::BufferUsageFlagBits::eStorageBuffer;
    buffer_info.sharingMode = vk::SharingMode::eExclusive;

    buffer_vk_ = vulkan_context.device().createBuffer(buffer_info);

    vk::MemoryRequirements mem_reqs = vulkan_context.device().getBufferMemoryRequirements(buffer_vk_);

    vk::MemoryAllocateInfo alloc_info{};
    alloc_info.allocationSize = mem_reqs.size;
    alloc_info.memoryTypeIndex = vulkan_context.memory_type_index();

    buffer_memory_vk_ = vulkan_context.device().allocateMemory(alloc_info);

    vulkan_context.device().bindBufferMemory(buffer_vk_, buffer_memory_vk_, 0);

    T* data = static_cast<T*>(vulkan_context.device().mapMemory(buffer_memory_vk_, 0, sizeof(T) * rows_ * cols_, vk::MemoryMapFlags()));
    std::copy(mat.data(), mat.data() + rows_ * cols_, data);
    vulkan_context.device().unmapMemory(buffer_memory_vk_);
  }

  ~matrix_vk() {
    if (buffer_vk_) {
      vulkan_context.device().destroyBuffer(buffer_vk_);
    }
    if (buffer_memory_vk_) {
      vulkan_context.device().freeMemory(buffer_memory_vk_);
    }
  }
};

}  // namespace math
}  // namespace stan

#endif
#endif
