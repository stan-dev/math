#ifndef STAN_MATH_OPENCL_OPENCL_CONTEXT_HPP
#define STAN_MATH_OPENCL_OPENCL_CONTEXT_HPP
#ifdef STAN_OPENCL

#define DEVICE_FILTER CL_DEVICE_TYPE_ALL
#ifndef OPENCL_DEVICE_ID
#error OPENCL_DEVICE_ID_NOT_SET
#endif
#ifndef OPENCL_PLATFORM_ID
#error OPENCL_PLATFORM_ID_NOT_SET
#endif

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>

#include <CL/opencl.hpp>
#include <tbb/concurrent_vector.h>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <cerrno>

/** \ingroup opencl
 *  \defgroup opencl_context_group OpenCL Context
 *  Initialization for OpenCL Context:
 *    1. create context
 *    2. Find OpenCL platforms and devices available
 *    3. set up command queue
 *    4. set architecture dependent kernel parameters
 */
namespace stan {
namespace math {

/** \ingroup opencl_context_group
 * The <code>opencl_context_base</code> class represents an OpenCL context
 * in the standard Meyers singleton design pattern.
 *
 * See the OpenCL specification glossary for a list of terms:
 * https://www.khronos.org/registry/OpenCL/specs/opencl-1.2.pdf.
 * The context includes the set of devices available on the host, command
 * queues, manages kernels.
 *
 * This is designed so there's only one instance running on the host.
 *
 * Some design decisions that may need to be addressed later:
 * - we are assuming a single OpenCL platform. We may want to run on multiple
 * platforms simultaneously
 * - we are assuming a single OpenCL device. We may want to run on multiple
 * devices simultaneously
 */
class opencl_context_base {
  friend class opencl_context;

 private:
  /** \ingroup opencl_context_group
   * Construct the opencl_context by initializing the
   * OpenCL context, devices, command queues, and kernel
   * groups.
   *
   * This constructor does the following:
   * 1. Gets the available platforms and selects the platform
   *  with given id.
   * 2. Gets the available devices and selects the device with given id.
   * 3. Creates the OpenCL context with the device.
   * 4. Creates the OpenCL command queue for the selected device.
   * 5. Sets OpenCL device dependent kernel parameters
   *
   * @param platform_id id of the OpenCL platform to use
   * @param device_id id of the OpenCL device to use
   * @throw std::system_error if an OpenCL error occurs.
   */
  opencl_context_base(int platform_id = OPENCL_PLATFORM_ID,
                      int device_id = OPENCL_DEVICE_ID) {
    try {
      // platform
      cl::Platform::get(&platforms_);
      if (platform_id >= platforms_.size()) {
        system_error("OpenCL Initialization", "[Platform]", -1,
                     "CL_INVALID_PLATFORM");
      }
      platform_.push_back(platforms_[platform_id]);
      platform_name_ = platform_[0].getInfo<CL_PLATFORM_NAME>();
      platform_[0].getDevices(DEVICE_FILTER, &devices_);
      if (devices_.size() == 0) {
        system_error("OpenCL Initialization", "[Device]", -1,
                     "CL_DEVICE_NOT_FOUND");
      }
      if (device_id >= devices_.size()) {
        system_error("OpenCL Initialization", "[Device]", -1,
                     "CL_INVALID_DEVICE");
      }
      device_.push_back(devices_[device_id]);
      // context and queue
      cl_command_queue_properties device_properties;
      device_[0].getInfo<cl_command_queue_properties>(
          CL_DEVICE_QUEUE_PROPERTIES, &device_properties);
      device_[0].getInfo<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE,
                                 &max_thread_block_size_);
      std::vector<size_t> max_wg_sizes
          = device_[0].getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
      if (max_wg_sizes.size() < 3) {
        system_error("OpenCL Initialization", "[Device]", -1,
                     "The device does not support 3D work groups!");
      }

      context_ = cl::Context(device_[0]);
      if (device_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE) {
        command_queue_
            = cl::CommandQueue(context_, device_[0],
                               CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, nullptr);
        in_order_ = CL_FALSE;
      } else {
        command_queue_ = cl::CommandQueue(context_, device_[0], 0, nullptr);
        in_order_ = CL_TRUE;
      }
      int max_square_block_size
          = std::min({max_wg_sizes[0], max_wg_sizes[1],
                      static_cast<size_t>(std::sqrt(max_thread_block_size_))});

      // Does a compile time check of the maximum allowed
      // dimension of a square thread block size
      // WG size of (32,32) works on all recent GPUs but would fail on some
      // older integrated GPUs or CPUs
      if (max_square_block_size < base_opts_["THREAD_BLOCK_SIZE"]) {
        base_opts_["THREAD_BLOCK_SIZE"] = max_square_block_size;
        base_opts_["WORK_PER_THREAD"] = 1;
      }
      if (std::min(max_thread_block_size_, max_wg_sizes[0])
          < base_opts_["LOCAL_SIZE_"]) {
        // must be a power of base_opts_["REDUCTION_STEP_SIZE"]
        const int p = std::log(max_thread_block_size_)
                      / std::log(base_opts_["REDUCTION_STEP_SIZE"]);
        base_opts_["LOCAL_SIZE_"]
            = std::pow(base_opts_["REDUCTION_STEP_SIZE"], p);
      }
      // Thread block size for the Cholesky
      // TODO(Steve): This should be tuned in a higher part of the stan language
      if (max_thread_block_size_ >= 256) {
        tuning_opts_.cholesky_min_L11_size = 256;
      } else {
        tuning_opts_.cholesky_min_L11_size = max_thread_block_size_;
      }
    } catch (const cl::Error& e) {
      check_opencl_error("opencl_context", e);
    }
  }

 protected:
  cl::Context context_;  // Manages the the device, queue, platform, memory, etc
  cl::CommandQueue command_queue_;       // job queue for device, one per device
  std::vector<cl::Platform> platforms_;  // Vector of available platforms
  std::vector<cl::Platform> platform_;   // The platform for compiling kernels
  std::string platform_name_;  // The platform such as NVIDIA OpenCL or AMD SDK
  std::vector<cl::Device> device_;   // The selected OpenCL device
  std::vector<cl::Device> devices_;  // All available OpenCL devices
  std::string device_name_;          // The name of OpenCL device
  size_t max_thread_block_size_;  // The maximum size of a block of workers on
                                  // the device
  bool in_order_;                 // Whether to use out of order execution.
  // Holds Default parameter values for each Kernel.
  using map_base_opts = std::unordered_map<std::string, int>;
  map_base_opts base_opts_
      = {{"LOWER", static_cast<int>(matrix_cl_view::Lower)},
         {"UPPER", static_cast<int>(matrix_cl_view::Upper)},
         {"ENTIRE", static_cast<int>(matrix_cl_view::Entire)},
         {"DIAGONAL", static_cast<int>(matrix_cl_view::Diagonal)},
         {"THREAD_BLOCK_SIZE", 32},
         {"WORK_PER_THREAD", 8},
         {"REDUCTION_STEP_SIZE", 4},
         {"LOCAL_SIZE_", 4096}};
  // TODO(Steve): Make these tunable during warmup
  struct tuning_struct {
    // Used in math/opencl/cholesky_decompose
    int cholesky_min_L11_size = 256;
    int cholesky_partition = 4;
    int cholesky_size_worth_transfer = 1250;
    // Used in math/rev/fun/cholesky_decompose
    int cholesky_rev_min_block_size = 512;
    int cholesky_rev_block_partition = 8;
    // used in math/opencl/multiply
    int multiply_wgs_per_compute_unit = 5;
    // used in math/prim/fun/gp_exp_quad_cov
    double gp_exp_quad_cov_complex = 1'000'000;
    double gp_exp_quad_cov_simple = 1'250;
    // used in math/prim/fun/multiply
    // and math/rev/fun/multiply
    int multiply_dim_prod_worth_transfer = 2000000;
    // used in math/prim/fun/mdivide_left_tri
    // and math/rev/fun/mdivide_left_tri
    int tri_inverse_size_worth_transfer = 100;
  } tuning_opts_;

 protected:
  static opencl_context_base& getInstance() noexcept {
    static opencl_context_base instance_;
    return instance_;
  }

  static void select_device(int platform_id, int device_id) {
    getInstance() = opencl_context_base(platform_id, device_id);
  }
};

/** \ingroup opencl_context_group
 * The API to access the methods and values in opencl_context_base
 */
class opencl_context {
  tbb::concurrent_vector<cl::Kernel*> kernel_caches_;

 public:
  opencl_context() = default;

  /** \ingroup opencl_context_group
   * Returns the description of the OpenCL platform and device that is used.
   * Devices will be an OpenCL and Platforms are a specific OpenCL
   * implementation such as AMD SDK's or Nvidia's OpenCL implementation.
   */
  inline std::string description() const {
    std::ostringstream msg;

    msg << "Platform ID: " << OPENCL_DEVICE_ID << "\n";
    msg << "Platform Name: "
        << opencl_context_base::getInstance()
               .platform_[0]
               .getInfo<CL_PLATFORM_NAME>()
        << "\n";
    msg << "Platform Vendor: "
        << opencl_context_base::getInstance()
               .platform_[0]
               .getInfo<CL_PLATFORM_VENDOR>()
        << "\n";
    msg << "\tDevice " << OPENCL_DEVICE_ID << ": "
        << "\n";
    msg << "\t\tDevice Name: "
        << opencl_context_base::getInstance()
               .device_[0]
               .getInfo<CL_DEVICE_NAME>()
        << "\n";
    msg << "\t\tDevice Type: "
        << opencl_context_base::getInstance()
               .device_[0]
               .getInfo<CL_DEVICE_TYPE>()
        << "\n";
    msg << "\t\tDevice Vendor: "
        << opencl_context_base::getInstance()
               .device_[0]
               .getInfo<CL_DEVICE_VENDOR>()
        << "\n";
    msg << "\t\tDevice Max Compute Units: "
        << opencl_context_base::getInstance()
               .device_[0]
               .getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()
        << "\n";
    msg << "\t\tDevice Global Memory: "
        << opencl_context_base::getInstance()
               .device_[0]
               .getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()
        << "\n";
    msg << "\t\tDevice Max Clock Frequency: "
        << opencl_context_base::getInstance()
               .device_[0]
               .getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>()
        << "\n";
    msg << "\t\tDevice Max Allocateable Memory: "
        << opencl_context_base::getInstance()
               .device_[0]
               .getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>()
        << "\n";
    msg << "\t\tDevice Local Memory: "
        << opencl_context_base::getInstance()
               .device_[0]
               .getInfo<CL_DEVICE_LOCAL_MEM_SIZE>()
        << "\n";
    msg << "\t\tDevice Available: "
        << opencl_context_base::getInstance()
               .device_[0]
               .getInfo<CL_DEVICE_AVAILABLE>()
        << "\n";
    return msg.str();
  }

  /** \ingroup opencl_context_group
   * Returns the description of the OpenCL platforms and devices that
   * are available. Devices will be an OpenCL and Platforms are a specific
   * OpenCL implementation such as AMD SDK's or Nvidia's OpenCL implementation.
   */
  inline std::string capabilities() const {
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);
    std::ostringstream msg;
    int platform_id = 0;
    int device_id = 0;

    msg << "Number of Platforms: " << all_platforms.size() << "\n";
    for (auto plat_iter : all_platforms) {
      cl::Platform platform(plat_iter);

      msg << "Platform ID: " << platform_id++ << "\n";
      msg << "Platform Name: " << platform.getInfo<CL_PLATFORM_NAME>() << "\n";
      msg << "Platform Vendor: " << platform.getInfo<CL_PLATFORM_VENDOR>()
          << "\n";

      try {
        std::vector<cl::Device> all_devices;
        platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);

        for (auto device_iter : all_devices) {
          cl::Device device(device_iter);

          msg << "\tDevice " << device_id++ << ": "
              << "\n";
          msg << "\t\tDevice Name: " << device.getInfo<CL_DEVICE_NAME>()
              << "\n";
          msg << "\t\tDevice Type: " << device.getInfo<CL_DEVICE_TYPE>()
              << "\n";
          msg << "\t\tDevice Vendor: " << device.getInfo<CL_DEVICE_VENDOR>()
              << "\n";
          msg << "\t\tDevice Max Compute Units: "
              << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << "\n";
          msg << "\t\tDevice Global Memory: "
              << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << "\n";
          msg << "\t\tDevice Max Clock Frequency: "
              << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << "\n";
          msg << "\t\tDevice Max Allocateable Memory: "
              << device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << "\n";
          msg << "\t\tDevice Local Memory: "
              << device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << "\n";
          msg << "\t\tDevice Available: "
              << device.getInfo<CL_DEVICE_AVAILABLE>() << "\n";
        }
      } catch (const cl::Error& e) {
        // if one of the platforms have no devices that match the device type
        // it will throw the error == -1 (DEVICE_NOT_FOUND)
        // other errors will throw a system error
        if (e.err() == -1) {
          msg << "\tno (OpenCL) devices in the platform with ID " << platform_id
              << "\n";
        } else {
          check_opencl_error("capabilities", e);
        }
      }
    }
    return msg.str();
  }

  /** \ingroup opencl_context_group
   * Returns the reference to the OpenCL context. The OpenCL context manages
   * objects such as the device, memory, command queue, program, and kernel
   * objects. For stan, there should only be one context, queue, device, and
   * program with multiple kernels.
   */
  inline cl::Context& context() noexcept {
    return opencl_context_base::getInstance().context_;
  }
  /** \ingroup opencl_context_group
   * Returns the reference to the active OpenCL command queue for the device.
   * One command queue will exist per device where
   * kernels are placed on the command queue and by default executed in order.
   */
  inline cl::CommandQueue& queue() noexcept {
    return opencl_context_base::getInstance().command_queue_;
  }
  /** \ingroup opencl_context_group
   * Returns a copy of the map of kernel defines
   */
  inline opencl_context_base::map_base_opts& base_opts() noexcept {
    return opencl_context_base::getInstance().base_opts_;
  }
  /** \ingroup opencl_context_group
   * Returns the maximum thread block size defined by
   * CL_DEVICE_MAX_WORK_GROUP_SIZE for the device in the context. This is the
   * maximum product of thread block dimensions for a particular device. IE a
   * max workgroup of 256 would allow thread blocks of sizes (16,16), (128,2),
   * (8, 32), etc.
   */
  inline int max_thread_block_size() noexcept {
    return opencl_context_base::getInstance().max_thread_block_size_;
  }

  /** \ingroup opencl_context_group
   * Returns the thread block size for the Cholesky Decompositions L_11.
   */
  inline opencl_context_base::tuning_struct& tuning_opts() noexcept {
    return opencl_context_base::getInstance().tuning_opts_;
  }

  /** \ingroup opencl_context_group
   * Returns a vector containing the OpenCL device used to create the context
   */
  inline std::vector<cl::Device>& device() noexcept {
    return opencl_context_base::getInstance().device_;
  }

  /** \ingroup opencl_context_group
   * Returns a vector containing the OpenCL platform used to create the context
   */
  inline std::vector<cl::Platform>& platform() noexcept {
    return opencl_context_base::getInstance().platform_;
  }
  /** \ingroup opencl_context_group
   * Return a bool representing whether the write to the OpenCL device are
   * blocking
   */
  inline bool& in_order() noexcept {
    return opencl_context_base::getInstance().in_order_;
  }

  /**
   * Selects the OpenCL device to use from now on.
   *
   * No `matrix_cl` objects or created before this call should be used after the
   * call (including any that might be on the AD stack)!
   * @param platform_id id of the platform the device is part of
   * @param instance_id if of the device
   */
  inline void select_device(int platform_id, int instance_id) {
    for (cl::Kernel* cache : kernel_caches_) {
      *cache = cl::Kernel();
    }
    kernel_caches_.clear();
    opencl_context_base::select_device(platform_id, instance_id);
  }

  /**
   * Registers a cached kernel. The cache will be invalidated if a new OpenCL
   * device is selected.
   * @param cache pointer to a cached kernel.
   */
  inline void register_kernel_cache(cl::Kernel* cache) {
    kernel_caches_.push_back(cache);
  }
};
static opencl_context opencl_context;
}  // namespace math
}  // namespace stan

#endif
#endif
