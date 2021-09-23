#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 200

#include <CL/opencl.hpp>
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

const int numElements = 32;

int main(void)
{
    // Filter for a 2.0 platform and set it as the default
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    cl::Platform plat;
    for (auto &p : platforms) {
        std::string platver = p.getInfo<CL_PLATFORM_VERSION>();
        if (platver.find("OpenCL 2.") != std::string::npos) {
            plat = p;
        }
    }
    if (plat() == 0) {
        std::cout << "No OpenCL 2.0 platform found.";
        return -1;
    }

    cl::Platform newP = cl::Platform::setDefault(plat);
    if (newP != plat) {
        std::cout << "Error setting default platform.";
        return -1;
    }

    // Use C++11 raw string literals for kernel source code
    std::string kernel1{R"CLC(
        global int globalA;
        kernel void updateGlobal()
        {
          globalA = 75;
        }
    )CLC"};
    std::string kernel2{R"CLC(
        typedef struct { global int *bar; } Foo;
        kernel void vectorAdd(global const Foo* aNum, global const int *inputA, global const int *inputB,
                              global int *output, int val, write_only pipe int outPipe, queue_t childQueue)
        {
          output[get_global_id(0)] = inputA[get_global_id(0)] + inputB[get_global_id(0)] + val + *(aNum->bar);
          write_pipe(outPipe, &val);
          queue_t default_queue = get_default_queue();
          ndrange_t ndrange = ndrange_1D(get_global_size(0)/2, get_global_size(0)/2);

          // Have a child kernel write into third quarter of output
          enqueue_kernel(default_queue, CLK_ENQUEUE_FLAGS_WAIT_KERNEL, ndrange,
            ^{
                output[get_global_size(0)*2 + get_global_id(0)] =
                  inputA[get_global_size(0)*2 + get_global_id(0)] + inputB[get_global_size(0)*2 + get_global_id(0)] + globalA;
            });

          // Have a child kernel write into last quarter of output
          enqueue_kernel(childQueue, CLK_ENQUEUE_FLAGS_WAIT_KERNEL, ndrange,
            ^{
                output[get_global_size(0)*3 + get_global_id(0)] =
                  inputA[get_global_size(0)*3 + get_global_id(0)] + inputB[get_global_size(0)*3 + get_global_id(0)] + globalA + 2;
            });
        }
    )CLC"};

    // New simpler string interface style
    std::vector<std::string> programStrings;
    programStrings.push_back(kernel1);
    programStrings.push_back(kernel2);

    cl::Program vectorAddProgram(programStrings);
    try {
        vectorAddProgram.build("-cl-std=CL2.0");
    }
    catch (...) {
        // Print build info for all devices
        cl_int buildErr = CL_SUCCESS;
        auto buildInfo = vectorAddProgram.getBuildInfo<CL_PROGRAM_BUILD_LOG>(&buildErr);
        for (auto &pair : buildInfo) {
            std::cerr << pair.second << std::endl << std::endl;
        }

        return 1;
    }

    typedef struct { int *bar; } Foo;

    // Get and run kernel that initializes the program-scope global
    // A test for kernels that take no arguments
    auto program2Kernel =
        cl::KernelFunctor<>(vectorAddProgram, "updateGlobal");
    program2Kernel(
        cl::EnqueueArgs(
        cl::NDRange(1)));

    //////////////////
    // SVM allocations

    auto anSVMInt = cl::allocate_svm<int, cl::SVMTraitCoarse<>>();
    *anSVMInt = 5;
    cl::SVMAllocator<Foo, cl::SVMTraitCoarse<cl::SVMTraitReadOnly<>>> svmAllocReadOnly;
    auto fooPointer = cl::allocate_pointer<Foo>(svmAllocReadOnly);
    fooPointer->bar = anSVMInt.get();
    cl::SVMAllocator<int, cl::SVMTraitCoarse<>> svmAlloc;
    std::vector<int, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> inputA(numElements, 1, svmAlloc);
    cl::coarse_svm_vector<int> inputB(numElements, 2, svmAlloc);

    //
    //////////////

    // Traditional cl_mem allocations
    std::vector<int> output(numElements, 0xdeadbeef);
    cl::Buffer outputBuffer(begin(output), end(output), false);
    cl::Pipe aPipe(sizeof(cl_int), numElements / 2);

    // Default command queue, also passed in as a parameter
    cl::DeviceCommandQueue defaultDeviceQueue = cl::DeviceCommandQueue::makeDefault(
        cl::Context::getDefault(), cl::Device::getDefault());

    auto vectorAddKernel =
        cl::KernelFunctor<
            decltype(fooPointer)&,
            int*,
            cl::coarse_svm_vector<int>&,
            cl::Buffer,
            int,
            cl::Pipe&,
            cl::DeviceCommandQueue
            >(vectorAddProgram, "vectorAdd");

    // Ensure that the additional SVM pointer is available to the kernel
    // This one was not passed as a parameter
    vectorAddKernel.setSVMPointers(anSVMInt);

	cl_int error;
	vectorAddKernel(
        cl::EnqueueArgs(
            cl::NDRange(numElements/2),
            cl::NDRange(numElements/2)),
        fooPointer,
        inputA.data(),
        inputB,
        outputBuffer,
        3,
        aPipe,
        defaultDeviceQueue,
		error
        );

    cl::copy(outputBuffer, begin(output), end(output));

    cl::Device d = cl::Device::getDefault();

    std::cout << "Output:\n";
    for (int i = 1; i < numElements; ++i) {
        std::cout << "\t" << output[i] << "\n";
    }
    std::cout << "\n\n";

    return 0;
}
