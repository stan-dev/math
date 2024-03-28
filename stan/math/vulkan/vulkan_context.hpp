#ifndef STAN_MATH_VULKAN_VULKAN_CONTEXT_HPP
#define STAN_MATH_VULKAN_VULKAN_CONTEXT_HPP
#ifdef STAN_VULKAN

#include <iostream>
#include <fstream>
#include <vulkan/vulkan.hpp>

namespace stan {
namespace math {

class vulkan_context_base {
  friend class vulkan_context;

  protected:
    vk::Instance instance_;
    vk::PhysicalDevice physical_device_;
    vk::PhysicalDeviceProperties device_properties_;
    uint32_t queue_family_index_ = 0;
    vk::Device device_;
    vk::PhysicalDeviceMemoryProperties memory_properties_;
    vk::CommandPool command_pool_;
    uint32_t memory_type_index_ = 0;
    vk::DeviceSize memory_heap_size_ = 0;
    vk::Queue queue_;
    size_t max_thread_block_size_;

    static vulkan_context_base& getInstance() noexcept {
      static vulkan_context_base context_base_;
      return context_base_;
    }

    static void select_device(int device_id) {
      getInstance() = vulkan_context_base(device_id);
    }

  private:
  vulkan_context_base(int device_id = VULKAN_PLATFORM_ID) {
    initialiseInstance();
    initialisePhysicalDevice(device_id);
    initialiseLogicalDevice();
  }

  void initialiseInstance() {
    vk::ApplicationInfo AppInfo{};
    AppInfo.pApplicationName = "Stan";
    AppInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    AppInfo.pEngineName = "Stan";
    AppInfo.engineVersion = VK_MAKE_VERSION(1, 0, 0);
    AppInfo.apiVersion = VK_API_VERSION_1_0;

    std::vector<const char*> requiredExtensions;
    vk::InstanceCreateInfo InstanceCreateInfo{};

  #ifdef __APPLE__
    requiredExtensions.emplace_back(VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME);
    InstanceCreateInfo.flags |= vk::InstanceCreateFlagBits::eEnumeratePortabilityKHR;
  #endif
    InstanceCreateInfo.pApplicationInfo = &AppInfo;
    InstanceCreateInfo.enabledLayerCount = 0;
    InstanceCreateInfo.enabledExtensionCount = requiredExtensions.size();
    InstanceCreateInfo.ppEnabledExtensionNames = requiredExtensions.data();

    instance_ = vk::createInstance(InstanceCreateInfo);
  }

  void initialisePhysicalDevice(int device_id) {
    std::vector<vk::PhysicalDevice> devices = instance_.enumeratePhysicalDevices();
    physical_device_ = devices[device_id];
    device_properties_ = physical_device_.getProperties();
    max_thread_block_size_ = device_properties_.limits.maxComputeWorkGroupInvocations;


    memory_properties_ = physical_device_.getMemoryProperties();

    for (uint32_t CurrentMemoryTypeIndex = 0; CurrentMemoryTypeIndex < memory_properties_.memoryTypeCount; ++CurrentMemoryTypeIndex) {
      vk::MemoryType MemoryType = memory_properties_.memoryTypes[CurrentMemoryTypeIndex];
      if ((vk::MemoryPropertyFlagBits::eHostVisible & MemoryType.propertyFlags) &&
        (vk::MemoryPropertyFlagBits::eHostCoherent & MemoryType.propertyFlags)) {
        memory_heap_size_ = memory_properties_.memoryHeaps[MemoryType.heapIndex].size;
        memory_type_index_ = CurrentMemoryTypeIndex;
        break;
      }
    }

    std::vector<vk::QueueFamilyProperties> QueueFamilyProps = physical_device_.getQueueFamilyProperties();
    auto PropIt = std::find_if(
      QueueFamilyProps.begin(), QueueFamilyProps.end(), [](const vk::QueueFamilyProperties& Prop) {
      return Prop.queueFlags & vk::QueueFlagBits::eCompute;
    });
    queue_family_index_ = std::distance(QueueFamilyProps.begin(), PropIt);
  }

  void initialiseLogicalDevice() {
    // Just to avoid a warning from the Vulkan Validation Layer
    const float QueuePriority = 1.0f;

    vk::DeviceQueueCreateInfo DeviceQueueCreateInfo{};
    DeviceQueueCreateInfo.flags |= vk::DeviceQueueCreateFlags();
    DeviceQueueCreateInfo.queueFamilyIndex = queue_family_index_;
    DeviceQueueCreateInfo.queueCount = 1;
    DeviceQueueCreateInfo.pQueuePriorities = &QueuePriority;

    vk::DeviceCreateInfo DeviceCreateInfo(vk::DeviceCreateFlags(),
                                          DeviceQueueCreateInfo);
    device_ = physical_device_.createDevice(DeviceCreateInfo);
    queue_ = device_.getQueue(queue_family_index_, 0);

    vk::CommandPoolCreateInfo CommandPoolCreateInfo(vk::CommandPoolCreateFlags(), queue_family_index_);
    command_pool_ = device_.createCommandPool(CommandPoolCreateInfo);
  }

  ~vulkan_context_base() {
    device_.destroyCommandPool(command_pool_);
    device_.destroy();
    instance_.destroy();
  }
};

class vulkan_context {
  public:
    vulkan_context() = default;

    inline vk::Instance& instance() noexcept {
      return vulkan_context_base::getInstance().instance_;
    }
    inline vk::PhysicalDevice& physical_device() noexcept {
      return vulkan_context_base::getInstance().physical_device_;
    }
    inline vk::Device& device() noexcept {
      return vulkan_context_base::getInstance().device_;
    }
    inline vk::CommandPool& command_pool() noexcept {
      return vulkan_context_base::getInstance().command_pool_;
    }
    inline vk::Queue& queue() noexcept {
      return vulkan_context_base::getInstance().queue_;
    }
    inline uint32_t& memory_type_index() noexcept {
      return vulkan_context_base::getInstance().memory_type_index_;
    }
    inline int max_thread_block_size() noexcept {
      return vulkan_context_base::getInstance().max_thread_block_size_;
    }

    void printProperties() noexcept {
      const auto device_properties_ = vulkan_context_base::getInstance().device_properties_;

      const uint32_t ApiVersion = device_properties_.apiVersion;
      vk::PhysicalDeviceLimits DeviceLimits = device_properties_.limits;

      std::cout
        << "Device Name: " << device_properties_.deviceName
        << "\nVulkan Version: " << VK_VERSION_MAJOR(ApiVersion) << "."
                                  << VK_VERSION_MINOR(ApiVersion) << "."
                                  << VK_VERSION_PATCH(ApiVersion)
        << "\nMax Compute Shared Memory Size: " << DeviceLimits.maxComputeSharedMemorySize / 1024 << " KB"
        << "\nCompute Queue Family Index: " << vulkan_context_base::getInstance().queue_family_index_
        << "\nMemory Type Index: " << vulkan_context_base::getInstance().memory_type_index_
        << "\nMemory Heap Size : " << vulkan_context_base::getInstance().memory_heap_size_ / 1024 / 1024 / 1024 << " GB"
        << std::endl;
    }
};

static vulkan_context vulkan_context;

}  // namespace math
}  // namespace stan

#endif
#endif
