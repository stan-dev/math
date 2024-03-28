#ifndef STAN_MATH_VULKAN_VULKAN_OPERATION_HPP
#define STAN_MATH_VULKAN_VULKAN_OPERATION_HPP
#ifdef STAN_VULKAN

#include <stan/math/vulkan/vulkan_context.hpp>
#include <stan/math/vulkan/matrix_vk.hpp>
#include <fstream>
#include <iostream>

namespace stan {
namespace math {

class vulkan_operation {
  vk::ShaderModule ShaderModule;
  vk::DescriptorSetLayout DescriptorSetLayout;
  vk::DescriptorPool DescriptorPool;
  vk::DescriptorSet DescriptorSet;
  vk::Pipeline ComputePipeline;
  vk::PipelineLayout PipelineLayout;
  vk::PipelineCache PipelineCache;
  vk::CommandBuffer CommandBuffer;

  private:
    void initialisePipeline() {
      vk::PipelineLayoutCreateInfo PipelineLayoutCreateInfo{};
      PipelineLayoutCreateInfo.setLayoutCount = 1;
      PipelineLayoutCreateInfo.pSetLayouts = &DescriptorSetLayout;

      PipelineLayout = vulkan_context.device().createPipelineLayout(PipelineLayoutCreateInfo);
      PipelineCache = vulkan_context.device().createPipelineCache(vk::PipelineCacheCreateInfo());

      vk::PipelineShaderStageCreateInfo ShaderStageCreateInfo{};
      ShaderStageCreateInfo.stage = vk::ShaderStageFlagBits::eCompute;
      ShaderStageCreateInfo.module = ShaderModule;
      ShaderStageCreateInfo.pName = "add";

      vk::ComputePipelineCreateInfo ComputePipelineCreateInfo{};
      ComputePipelineCreateInfo.stage = ShaderStageCreateInfo;
      ComputePipelineCreateInfo.layout = PipelineLayout;

      ComputePipeline = vulkan_context.device().createComputePipeline(PipelineCache, ComputePipelineCreateInfo).value;
    }

    void initialiseCommandBuffer() {
      vk::CommandBufferAllocateInfo CommandBufferAllocInfo{};
      CommandBufferAllocInfo.commandPool = vulkan_context.command_pool();
      CommandBufferAllocInfo.level = vk::CommandBufferLevel::ePrimary;
      CommandBufferAllocInfo.commandBufferCount = 1;

      CommandBuffer = vulkan_context.device().allocateCommandBuffers(CommandBufferAllocInfo).front();

      vk::CommandBufferBeginInfo CommandBufferBeginInfo{};
      CommandBuffer.begin(CommandBufferBeginInfo);

      CommandBuffer.bindPipeline(vk::PipelineBindPoint::eCompute, ComputePipeline);
      CommandBuffer.bindDescriptorSets(vk::PipelineBindPoint::eCompute, PipelineLayout, 0, DescriptorSet, {});
      CommandBuffer.dispatch(vulkan_context.max_thread_block_size(), 1, 1);
      CommandBuffer.end();
    }

    void submitCommandBuffer() {
      vk::SubmitInfo SubmitInfo{};
      SubmitInfo.commandBufferCount = 1;
      SubmitInfo.pCommandBuffers = &CommandBuffer;

      vk::Fence Fence = vulkan_context.device().createFence(vk::FenceCreateInfo());
      vulkan_context.queue().submit(SubmitInfo, Fence);
      vulkan_context.device().waitForFences(Fence, VK_TRUE, std::numeric_limits<uint64_t>::max());
    }

  public:
    vulkan_operation(const std::string& filename) {
      std::vector<char> ShaderContents;
      if (std::ifstream ShaderFile{ filename, std::ios::binary | std::ios::ate }) {
        const size_t FileSize = ShaderFile.tellg();
        ShaderFile.seekg(0);
        ShaderContents.resize(FileSize, '\0');
        ShaderFile.read(ShaderContents.data(), FileSize);
      }
      vk::ShaderModuleCreateInfo ShaderModuleCreateInfo{};
      ShaderModuleCreateInfo.flags = vk::ShaderModuleCreateFlags();
      ShaderModuleCreateInfo.codeSize = ShaderContents.size();
      ShaderModuleCreateInfo.pCode = reinterpret_cast<const uint32_t*>(ShaderContents.data());

      ShaderModule = vulkan_context.device().createShaderModule(ShaderModuleCreateInfo);
    }

    ~vulkan_operation() {
      vulkan_context.device().destroyPipeline(ComputePipeline);
      vulkan_context.device().destroyPipelineLayout(PipelineLayout);
      vulkan_context.device().destroyPipelineCache(PipelineCache);
      vulkan_context.device().destroyShaderModule(ShaderModule);
      vulkan_context.device().destroyDescriptorSetLayout(DescriptorSetLayout);
      vulkan_context.device().destroyDescriptorPool(DescriptorPool);
    }

    template<typename MatrixVkResultT, typename... MatrixVkInputsT>
    void add_operands(MatrixVkResultT&& matrix_vk_result, MatrixVkInputsT&&... matrix_vk_inputs) {
      std::vector<vk::DescriptorSetLayoutBinding> DescriptorSetLayoutBindings;
      std::vector<vk::DescriptorBufferInfo> DescriptorBufferInfos;
      std::vector<vk::WriteDescriptorSet> WriteDescriptorSets;

      int Binding = 0;
      auto add_operand = [&](auto& matrix_vk) {
        vk::DescriptorSetLayoutBinding DescriptorSetLayoutBinding{};
        DescriptorSetLayoutBinding.binding = Binding;
        DescriptorSetLayoutBinding.descriptorType = vk::DescriptorType::eStorageBuffer;
        DescriptorSetLayoutBinding.descriptorCount = 1;
        DescriptorSetLayoutBinding.stageFlags = vk::ShaderStageFlagBits::eCompute;
        DescriptorSetLayoutBinding.pImmutableSamplers = nullptr;
        DescriptorSetLayoutBindings.push_back(DescriptorSetLayoutBinding);

        vk::DescriptorBufferInfo DescriptorBufferInfo{};
        DescriptorBufferInfo.buffer = matrix_vk.buffer();
        DescriptorBufferInfo.offset = 0;
        DescriptorBufferInfo.range = VK_WHOLE_SIZE;
        DescriptorBufferInfos.push_back(DescriptorBufferInfo);

        ++Binding;
        return 0;
      };

      add_operand(matrix_vk_result);
      static_cast<void>(std::initializer_list<int>{ 0, add_operand(matrix_vk_inputs)... });


      vk::DescriptorSetLayoutCreateInfo DescriptorSetLayoutCreateInfo(vk::DescriptorSetLayoutCreateFlags(),
                                      DescriptorSetLayoutBindings);
      DescriptorSetLayout = vulkan_context.device().createDescriptorSetLayout(DescriptorSetLayoutCreateInfo);

      vk::DescriptorPoolSize DescriptorPoolSize{};
      DescriptorPoolSize.type = vk::DescriptorType::eStorageBuffer;
      DescriptorPoolSize.descriptorCount = DescriptorBufferInfos.size();

      vk::DescriptorPoolCreateInfo DescriptorPoolCreateInfo{};
      DescriptorPoolCreateInfo.flags = vk::DescriptorPoolCreateFlags();
      DescriptorPoolCreateInfo.maxSets = 1;
      DescriptorPoolCreateInfo.poolSizeCount = 1;
      DescriptorPoolCreateInfo.pPoolSizes = &DescriptorPoolSize;

      DescriptorPool = vulkan_context.device().createDescriptorPool(DescriptorPoolCreateInfo);

      vk::DescriptorSetAllocateInfo DescriptorSetAllocInfo(DescriptorPool, 1, &DescriptorSetLayout);
      const std::vector<vk::DescriptorSet> DescriptorSets = vulkan_context.device().allocateDescriptorSets(DescriptorSetAllocInfo);
      DescriptorSet = DescriptorSets.front();

      Binding = 0;
      auto add_operand_write = [&](auto& matrix_vk) {
        vk::WriteDescriptorSet WriteDescriptorSet{};
        WriteDescriptorSet.dstSet = DescriptorSet;
        WriteDescriptorSet.dstBinding = Binding;
        WriteDescriptorSet.dstArrayElement = 0;
        WriteDescriptorSet.descriptorType = vk::DescriptorType::eStorageBuffer;
        WriteDescriptorSet.descriptorCount = 1;
        WriteDescriptorSet.pBufferInfo = &DescriptorBufferInfos[Binding];
        WriteDescriptorSets.push_back(WriteDescriptorSet);

        ++Binding;
        return 0;
      };

      add_operand_write(matrix_vk_result);
      static_cast<void>(std::initializer_list<int>{ 0, add_operand_write(matrix_vk_inputs)... });

      vulkan_context.device().updateDescriptorSets(WriteDescriptorSets, {});
    }

    void eval_operation() {
      initialisePipeline();
      initialiseCommandBuffer();
      submitCommandBuffer();
    }
};

}  // namespace math
}  // namespace stan
#endif
#endif
