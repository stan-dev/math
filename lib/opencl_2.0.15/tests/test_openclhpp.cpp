#define CL_HPP_UNIT_TEST_ENABLE
#define CL_HPP_USE_CL_SUB_GROUPS_KHR

// We want to support all versions
#define CL_HPP_MINIMUM_OPENCL_VERSION 100
# include <CL/opencl.hpp>
# define TEST_RVALUE_REFERENCES
# define VECTOR_CLASS cl::vector
# define STRING_CLASS cl::string

extern "C"
{
#include <unity.h>
#include <cmock.h>
#include "Mockcl.h"
#include <string.h>

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

/// Creates fake IDs that are easy to identify

static inline cl_platform_id make_platform_id(int index)
{
    return (cl_platform_id) (size_t) (0x1a1a1a1a + index);
}

static inline cl_context make_context(int index)
{
    return (cl_context) (size_t) (0xcccccccc + index);
}

static inline cl_device_id make_device_id(int index)
{
    return (cl_device_id) (size_t) (0xdededede + index);
}

static inline cl_mem make_mem(int index)
{
    return (cl_mem) (size_t) (0x33333333 + index);
}

static inline cl_command_queue make_command_queue(int index)
{
    return (cl_command_queue) (size_t) (0xc0c0c0c0 + index);
}

static inline cl_kernel make_kernel(int index)
{
    return (cl_kernel) (size_t) (0xcececece + index);
}

static inline cl_program make_program(int index)
{
    return (cl_program)(size_t)(0xcfcfcfcf + index);
}

/* Pools of pre-allocated wrapped objects for tests. There is no device pool,
 * because there is no way to know whether the test wants the device to be
 * reference countable or not.
 */
static const int POOL_MAX = 5;
static cl::Platform platformPool[POOL_MAX];
static cl::Context contextPool[POOL_MAX];
static cl::CommandQueue commandQueuePool[POOL_MAX];
static cl::Buffer bufferPool[POOL_MAX];
static cl::Image2D image2DPool[POOL_MAX];
static cl::Image3D image3DPool[POOL_MAX];
static cl::Kernel kernelPool[POOL_MAX];

/****************************************************************************
 * Stub functions shared by multiple tests
 ****************************************************************************/

/**
 * Stub implementation of clGetCommandQueueInfo that returns the first context.
 */
static cl_int clGetCommandQueueInfo_context(
    cl_command_queue id,
    cl_command_queue_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void)num_calls;


    TEST_ASSERT_EQUAL_HEX(CL_QUEUE_CONTEXT, param_name);
    TEST_ASSERT(param_value == NULL || param_value_size >= sizeof(cl_context));
    if (param_value_size_ret != NULL)
        *param_value_size_ret = sizeof(cl_context);
    if (param_value != NULL)
        *(cl_context *)param_value = make_context(0);

    return CL_SUCCESS;
}

/**
 * Stub implementation of clGetDeviceInfo that just returns the first platform.
 */
static cl_int clGetDeviceInfo_platform(
    cl_device_id id,
    cl_device_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void) num_calls;

    TEST_ASSERT_EQUAL_HEX(CL_DEVICE_PLATFORM, param_name);
    TEST_ASSERT(param_value == NULL || param_value_size >= sizeof(cl_platform_id));
    if (param_value_size_ret != NULL)
        *param_value_size_ret = sizeof(cl_platform_id);
    if (param_value != NULL)
        *(cl_platform_id *) param_value = make_platform_id(0);
    return CL_SUCCESS;
}

/**
 * Stub implementation of clGetContextInfo that just returns the first device.
 */
static cl_int clGetContextInfo_device(
    cl_context id,
    cl_context_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void) num_calls;

    TEST_ASSERT_EQUAL_HEX(CL_CONTEXT_DEVICES, param_name);
    TEST_ASSERT(param_value == NULL || param_value_size >= sizeof(cl_device_id));
    if (param_value_size_ret != NULL)
        *param_value_size_ret = sizeof(cl_device_id);
    if (param_value != NULL)
        *(cl_device_id *) param_value = make_device_id(0);
    return CL_SUCCESS;
}


/**
 * Stub implementation of clGetPlatformInfo that returns a specific version.
 * It also checks that the id is the zeroth platform.
 */
static cl_int clGetPlatformInfo_version(
    cl_platform_id id,
    cl_platform_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    const char *version)
{
    size_t bytes = strlen(version) + 1;

    TEST_ASSERT_NOT_NULL(id);
    TEST_ASSERT_EQUAL_PTR(make_platform_id(0), id);
    TEST_ASSERT_EQUAL_HEX(CL_PLATFORM_VERSION, param_name);
    TEST_ASSERT(param_value == NULL || param_value_size >= bytes);
    if (param_value_size_ret != NULL)
        *param_value_size_ret = bytes;
    if (param_value != NULL)
        strcpy((char *) param_value, version);
    return CL_SUCCESS;
}

/**
 * A stub for clGetPlatformInfo that will only support querying
 * CL_PLATFORM_VERSION, and will return version 1.1.
 */
static cl_int clGetPlatformInfo_version_1_1(
    cl_platform_id id,
    cl_platform_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void) num_calls;
    return clGetPlatformInfo_version(
        id, param_name, param_value_size, param_value,
        param_value_size_ret, "OpenCL 1.1 Mock");
}

/**
 * A stub for clGetPlatformInfo that will only support querying
 * CL_PLATFORM_VERSION, and will return version 1.2.
 */
static cl_int clGetPlatformInfo_version_1_2(
    cl_platform_id id,
    cl_platform_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void) num_calls;
    return clGetPlatformInfo_version(
        id, param_name, param_value_size, param_value,
        param_value_size_ret, "OpenCL 1.2 Mock");
}

/**
 * A stub for clGetPlatformInfo that will only support querying
 * CL_PLATFORM_VERSION, and will return version 2.0.
 */
static cl_int clGetPlatformInfo_version_2_0(
    cl_platform_id id,
    cl_platform_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void) num_calls;
    return clGetPlatformInfo_version(
        id, param_name, param_value_size, param_value,
        param_value_size_ret, "OpenCL 2.0 Mock");
}

/**
 * A stub for clGetPlatformInfo that will only support querying
 * CL_PLATFORM_VERSION, and will return version 3.0.
 */
static cl_int clGetPlatformInfo_version_3_0(
    cl_platform_id id,
    cl_platform_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void) num_calls;
    return clGetPlatformInfo_version(
        id, param_name, param_value_size, param_value,
        param_value_size_ret, "OpenCL 3.0 Mock");
}

/* Simulated reference counts. The table points to memory held by the caller.
 * This makes things simpler in the common case of only one object to be
 * reference counted.
 */
class RefcountTable
{
private:
    int n; // number of objects
    void * const *objects; // object IDs
    int *refcounts;        // current refcounts

    int find(void *object)
    {
        int idx = 0;
        while (idx < n && objects[idx] != object)
            idx++;
        TEST_ASSERT(idx < n);
        TEST_ASSERT(refcounts[idx] > 0); // otherwise object has been destroyed
        return idx;
    }

public:
    RefcountTable() : n(0), objects(NULL), refcounts(NULL) {}

    void init(int n, void * const *objects, int *refcounts)
    {
        this->n = n;
        this->objects = objects;
        this->refcounts = refcounts;
    }

    void reset()
    {
        init(0, NULL, NULL);
    }

    cl_int retain(void *object)
    {
        int idx = find(object);
        ++refcounts[idx];
        return CL_SUCCESS;
    }

    cl_int release(void *object)
    {
        int idx = find(object);
        --refcounts[idx];
        return CL_SUCCESS;
    }
};

/* Stubs for retain/release calls that track reference counts. The stubs
 * check that the reference count never becomes negative and that a zero
 * reference count is never incremented.
 *
 * Use the prepareRefcount* calls to set up the global variables first.
 */

#define MAKE_REFCOUNT_STUBS(cl_type, retainfunc, releasefunc, table) \
    static RefcountTable table; \
    static cl_int retainfunc ## _refcount(cl_type object, int num_calls) \
    { \
        (void) num_calls; \
        return table.retain(object); \
    } \
    static cl_int releasefunc ## _refcount(cl_type object, int num_calls) \
    { \
        (void) num_calls; \
        return table.release(object); \
    } \
    static void prepare_ ## table(int n, cl_type const *objects, int *refcounts) \
    { \
        table.init(n, (void * const *) objects, refcounts); \
        retainfunc ## _StubWithCallback(retainfunc ## _refcount); \
        releasefunc ## _StubWithCallback(releasefunc ## _refcount); \
    }

MAKE_REFCOUNT_STUBS(cl_device_id, clRetainDevice, clReleaseDevice, deviceRefcounts)
MAKE_REFCOUNT_STUBS(cl_context, clRetainContext, clReleaseContext, contextRefcounts)
MAKE_REFCOUNT_STUBS(cl_mem, clRetainMemObject, clReleaseMemObject, memRefcounts)

/* The indirection through MAKE_MOVE_TESTS2 with a prefix parameter is to
 * prevent the simple-minded parser from Unity from identifying tests from the
 * macro value.
 */
#ifdef TEST_RVALUE_REFERENCES
#define MAKE_MOVE_TESTS2(prefix, type, makeFunc, releaseFunc, pool) \
    void prefix ## MoveAssign ## type ## NonNull() \
    { \
        releaseFunc ## _ExpectAndReturn(makeFunc(0), CL_SUCCESS); \
        pool[0] = std::move(pool[1]); \
        TEST_ASSERT_EQUAL_PTR(makeFunc(1), pool[0]()); \
        TEST_ASSERT_NULL(pool[1]()); \
    } \
    \
    void prefix ## MoveAssign ## type ## Null() \
    { \
        pool[0]() = NULL; \
        pool[0] = std::move(pool[1]); \
        TEST_ASSERT_EQUAL_PTR(makeFunc(1), pool[0]()); \
        TEST_ASSERT_NULL(pool[1]()); \
    } \
    \
    void prefix ## MoveConstruct ## type ## NonNull() \
    { \
        cl::type tmp(std::move(pool[0])); \
        TEST_ASSERT_EQUAL_PTR(makeFunc(0), tmp()); \
        TEST_ASSERT_NULL(pool[0]()); \
        tmp() = NULL; \
    } \
    \
    void prefix ## MoveConstruct ## type ## Null() \
    { \
        cl::type empty; \
        cl::type tmp(std::move(empty)); \
        TEST_ASSERT_NULL(tmp()); \
        TEST_ASSERT_NULL(empty()); \
    }
#else
#define MAKE_MOVE_TESTS2(prefix, type, makeFunc, releaseFunc, pool) \
    void prefix ## MoveAssign ## type ## NonNull() {} \
    void prefix ## MoveAssign ## type ## Null() {} \
    void prefix ## MoveConstruct ## type ## NonNull() {} \
    void prefix ## MoveConstruct ## type ## Null() {}
#endif // !TEST_RVALUE_REFERENCES
#define MAKE_MOVE_TESTS(type, makeFunc, releaseFunc, pool) \
    MAKE_MOVE_TESTS2(test, type, makeFunc, releaseFunc, pool)

void setUp()
{
    /* We reach directly into the objects rather than using assignment to
     * avoid the reference counting functions from being called.
     */
    for (int i = 0; i < POOL_MAX; i++)
    {
        platformPool[i]() = make_platform_id(i);
        contextPool[i]() = make_context(i);
        commandQueuePool[i]() = make_command_queue(i);
        bufferPool[i]() = make_mem(i);
        image2DPool[i]() = make_mem(i);
        image3DPool[i]() = make_mem(i);
        kernelPool[i]() = make_kernel(i);
    }

    deviceRefcounts.reset();
    contextRefcounts.reset();
    memRefcounts.reset();
}

void tearDown()
{
    /* Wipe out the internal state to avoid a release call being made */
    for (int i = 0; i < POOL_MAX; i++)
    {
        platformPool[i]() = NULL;
        contextPool[i]() = NULL;
        commandQueuePool[i]() = NULL;
        bufferPool[i]() = NULL;
        image2DPool[i]() = NULL;
        image3DPool[i]() = NULL;
        kernelPool[i]() = NULL;
    }
}

/****************************************************************************
 * Tests for cl::Context
 ****************************************************************************/

void testCopyContextNonNull()
{
    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);
    clRetainContext_ExpectAndReturn(make_context(1), CL_SUCCESS);

    contextPool[0] = contextPool[1];
    TEST_ASSERT_EQUAL_PTR(make_context(1), contextPool[0]());
}

void testMoveAssignContextNonNull();
void testMoveAssignContextNull();
void testMoveConstructContextNonNull();
void testMoveConstructContextNull();
MAKE_MOVE_TESTS(Context, make_context, clReleaseContext, contextPool)

/// Stub for querying CL_CONTEXT_DEVICES that returns two devices
static cl_int clGetContextInfo_testContextGetDevices(
    cl_context context,
    cl_context_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void) num_calls;
    TEST_ASSERT_EQUAL_PTR(make_context(0), context);
    TEST_ASSERT_EQUAL_HEX(CL_CONTEXT_DEVICES, param_name);
    TEST_ASSERT(param_value == NULL || param_value_size >= 2 * sizeof(cl_device_id));
    if (param_value_size_ret != NULL)
        *param_value_size_ret = 2 * sizeof(cl_device_id);
    if (param_value != NULL)
    {
        cl_device_id *devices = (cl_device_id *) param_value;
        devices[0] = make_device_id(0);
        devices[1] = make_device_id(1);
    }
    return CL_SUCCESS;
}

/// Test that queried devices are not refcounted
void testContextGetDevices1_1()
{
    clGetContextInfo_StubWithCallback(clGetContextInfo_testContextGetDevices);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);

    VECTOR_CLASS<cl::Device> devices = contextPool[0].getInfo<CL_CONTEXT_DEVICES>();
    TEST_ASSERT_EQUAL(2, devices.size());
    TEST_ASSERT_EQUAL_PTR(make_device_id(0), devices[0]());
    TEST_ASSERT_EQUAL_PTR(make_device_id(1), devices[1]());
}

/// Test that queried devices are correctly refcounted
void testContextGetDevices1_2()
{
    clGetContextInfo_StubWithCallback(clGetContextInfo_testContextGetDevices);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);

    clRetainDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);
    clRetainDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);

    VECTOR_CLASS<cl::Device> devices = contextPool[0].getInfo<CL_CONTEXT_DEVICES>();
    TEST_ASSERT_EQUAL(2, devices.size());
    TEST_ASSERT_EQUAL_PTR(make_device_id(0), devices[0]());
    TEST_ASSERT_EQUAL_PTR(make_device_id(1), devices[1]());

    // Prevent release in the destructor
    devices[0]() = NULL;
    devices[1]() = NULL;
}

// This is used to get a list of all platforms, so expect two calls
// First, return to say we have two platforms
// Then return the two platform id_s
static cl_int clGetPlatformIDs_testContextFromType(
    cl_uint num_entries,
    cl_platform_id *platforms,
    cl_uint *num_platforms,
    int num_calls)
{
    if (num_calls == 0)
    {
        TEST_ASSERT_NULL(platforms);
        TEST_ASSERT_NOT_NULL(num_platforms);
        *num_platforms = 2;
        return CL_SUCCESS;
    }
    else if (num_calls == 1)
    {
        TEST_ASSERT_NOT_NULL(platforms);
        TEST_ASSERT_EQUAL(2, num_entries);
        platforms[0] = make_platform_id(0);
        platforms[1] = make_platform_id(1);
        return CL_SUCCESS;
    }
    else
    {
        TEST_FAIL_MESSAGE("clGetPlatformIDs called too many times");
        return CL_INVALID_VALUE;
    }
}

// Expect three calls to this
// 1. Platform 1, we have no GPUs
// 2. Platform 2, we have two GPUs
// 3. Here are the two cl_device_id's
static cl_int clGetDeviceIDs_testContextFromType(
    cl_platform_id  platform,
    cl_device_type  device_type,
    cl_uint  num_entries,
    cl_device_id  *devices,
    cl_uint  *num_devices,
    int num_calls)
{
    if (num_calls == 0)
    {
        TEST_ASSERT_EQUAL_PTR(make_platform_id(0), platform);
        TEST_ASSERT_EQUAL(CL_DEVICE_TYPE_GPU, device_type);
        TEST_ASSERT_NOT_NULL(num_devices);
        return CL_DEVICE_NOT_FOUND;
    }
    else if (num_calls == 1)
    {
        TEST_ASSERT_EQUAL_PTR(make_platform_id(1), platform);
        TEST_ASSERT_EQUAL(CL_DEVICE_TYPE_GPU, device_type);
        TEST_ASSERT_NOT_NULL(num_devices);
        *num_devices = 2;
        return CL_SUCCESS;
    }
    else if (num_calls == 2)
    {
        TEST_ASSERT_EQUAL_PTR(make_platform_id(1), platform);
        TEST_ASSERT_EQUAL(CL_DEVICE_TYPE_GPU, device_type);
        TEST_ASSERT_EQUAL(2, num_entries);
        TEST_ASSERT_NOT_NULL(devices);
        devices[0] = make_device_id(0);
        devices[1] = make_device_id(1);
        return CL_SUCCESS;
    }
    else
    {
        TEST_FAIL_MESSAGE("clGetDeviceIDs called too many times");
        return CL_INVALID_VALUE;
    }
}

// Stub for clCreateContextFromType
// - expect platform 1 with GPUs and non-null properties
static cl_context clCreateContextFromType_testContextFromType(
    const cl_context_properties  *properties,
    cl_device_type  device_type,
    void  (CL_CALLBACK *pfn_notify) (const char *errinfo,
    const void  *private_info,
    size_t  cb,
    void  *user_data),
    void  *user_data,
    cl_int  *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(CL_DEVICE_TYPE_GPU, device_type);
#if !defined(__APPLE__) && !defined(__MACOS)
    TEST_ASSERT_NOT_NULL(properties);
    TEST_ASSERT_EQUAL(CL_CONTEXT_PLATFORM, properties[0]);
    TEST_ASSERT_EQUAL(make_platform_id(1), properties[1]);
#endif
    return make_context(0);
}

void testContextFromType()
{
#if !defined(__APPLE__) && !defined(__MACOS)
    clGetPlatformIDs_StubWithCallback(clGetPlatformIDs_testContextFromType);
    clGetDeviceIDs_StubWithCallback(clGetDeviceIDs_testContextFromType);

    // The opencl.hpp header will perform an extra retain here to be consistent
    // with other APIs retaining runtime-owned objects before releasing them
    clRetainDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);
    clRetainDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);

    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);

    // End of scope of vector of devices within constructor
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);
    clReleaseDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);
#endif

    clCreateContextFromType_StubWithCallback(clCreateContextFromType_testContextFromType);

    cl::Context context(CL_DEVICE_TYPE_GPU);
    TEST_ASSERT_EQUAL_PTR(make_context(0), context());

    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);
}

void testContextFromTypeNonNullProperties()
{
    clCreateContextFromType_StubWithCallback(clCreateContextFromType_testContextFromType);

    const cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)make_platform_id(1), 0 };
    cl::Context context(CL_DEVICE_TYPE_GPU, props);
    TEST_ASSERT_EQUAL_PTR(make_context(0), context());

    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);
}

static cl_context clCreateContext_testContextNonNullProperties(
    const cl_context_properties* properties,
    cl_uint num_devices,
    const cl_device_id* devices,
    void  (CL_CALLBACK *pfn_notify) (const char *errinfo, const void  *private_info, size_t  cb, void  *user_data),
    void  *user_data,
    cl_int  *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_NOT_NULL(properties);
    TEST_ASSERT_GREATER_THAN(0, num_devices);
    for (int i = 0; i < num_devices; i++) {
        TEST_ASSERT_EQUAL(make_device_id(i), devices[i]);
    }
    return make_context(0);
}

void testContextWithDeviceNonNullProperties()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);
    clCreateContext_StubWithCallback(clCreateContext_testContextNonNullProperties);
    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);

    const cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)make_platform_id(0), 0 };
    cl::Device device = cl::Device(make_device_id(0));

    cl::Context context(device, props);
    TEST_ASSERT_EQUAL_PTR(make_context(0), context());
}

void testContextWithDevicesNonNullProperties()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);
    clCreateContext_StubWithCallback(clCreateContext_testContextNonNullProperties);
    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);

    const cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)make_platform_id(0), 0 };
    cl::Device device0 = cl::Device(make_device_id(0));
    cl::Device device1 = cl::Device(make_device_id(1));

    cl::Context context({device0, device1}, props);
    TEST_ASSERT_EQUAL_PTR(make_context(0), context());
}

/****************************************************************************
 * Tests for cl::CommandQueue
 ****************************************************************************/

void testMoveAssignCommandQueueNonNull();
void testMoveAssignCommandQueueNull();
void testMoveConstructCommandQueueNonNull();
void testMoveConstructCommandQueueNull();
MAKE_MOVE_TESTS(CommandQueue, make_command_queue, clReleaseCommandQueue, commandQueuePool);

// Stub for clGetCommandQueueInfo that returns context 0
static cl_int clGetCommandQueueInfo_testCommandQueueGetContext(
    cl_command_queue command_queue,
    cl_command_queue_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void) num_calls;
    TEST_ASSERT_EQUAL_PTR(make_command_queue(0), command_queue);
    TEST_ASSERT_EQUAL_HEX(CL_QUEUE_CONTEXT, param_name);
    TEST_ASSERT(param_value == NULL || param_value_size >= sizeof(cl_context));
    if (param_value_size_ret != NULL)
        *param_value_size_ret = sizeof(cl_context);
    if (param_value != NULL)
        *(cl_context *) param_value = make_context(0);
    return CL_SUCCESS;
}

void testCommandQueueGetContext()
{
    cl_context expected = make_context(0);
    int refcount = 1;

    clGetCommandQueueInfo_StubWithCallback(clGetCommandQueueInfo_testCommandQueueGetContext);
    prepare_contextRefcounts(1, &expected, &refcount);

    cl::Context ctx = commandQueuePool[0].getInfo<CL_QUEUE_CONTEXT>();
    TEST_ASSERT_EQUAL_PTR(expected, ctx());
    TEST_ASSERT_EQUAL(2, refcount);

    ctx() = NULL;
}

// Stub for clGetCommandQueueInfo that returns device 0
static cl_int clGetCommandQueueInfo_testCommandQueueGetDevice(
    cl_command_queue command_queue,
    cl_command_queue_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void) num_calls;
    TEST_ASSERT_EQUAL_PTR(make_command_queue(0), command_queue);
    TEST_ASSERT_EQUAL_HEX(CL_QUEUE_DEVICE, param_name);
    TEST_ASSERT(param_value == NULL || param_value_size >= sizeof(cl_device_id));
    if (param_value_size_ret != NULL)
        *param_value_size_ret = sizeof(cl_device_id);
    if (param_value != NULL)
        *(cl_device_id *) param_value = make_device_id(0);
    return CL_SUCCESS;
}

void testCommandQueueGetDevice1_1()
{
    cl_device_id expected = make_device_id(0);

    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);
    clGetCommandQueueInfo_StubWithCallback(clGetCommandQueueInfo_testCommandQueueGetDevice);

    cl::Device device = commandQueuePool[0].getInfo<CL_QUEUE_DEVICE>();
    TEST_ASSERT_EQUAL_PTR(expected, device());

    device() = NULL;
}

void testCommandQueueGetDevice1_2()
{
    cl_device_id expected = make_device_id(0);
    int refcount = 1;

    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clGetCommandQueueInfo_StubWithCallback(clGetCommandQueueInfo_testCommandQueueGetDevice);
    prepare_deviceRefcounts(1, &expected, &refcount);

    cl::Device device = commandQueuePool[0].getInfo<CL_QUEUE_DEVICE>();
    TEST_ASSERT_EQUAL_PTR(expected, device());
    TEST_ASSERT_EQUAL(2, refcount);

    device() = NULL;
}

// stub for clCreateCommandQueue - returns queue zero
static cl_command_queue clCreateCommandQueue_testCommandQueueFromSpecifiedContext(
    cl_context context,
    cl_device_id device,
    cl_command_queue_properties properties,
    cl_int *errcode_ret,
    int num_calls)
{
    (void) num_calls;
    TEST_ASSERT_EQUAL_PTR(make_context(0), context);
    TEST_ASSERT_EQUAL_PTR(make_device_id(0), device);
    TEST_ASSERT(properties == 0);
    return make_command_queue(0);
}

#if CL_HPP_TARGET_OPENCL_VERSION >= 200
// stub for clCreateCommandQueueWithProperties - returns queue zero
static cl_command_queue clCreateCommandQueueWithProperties_testCommandQueueFromSpecifiedContext(
    cl_context context,
    cl_device_id device,
    const cl_queue_properties *properties,
    cl_int *errcode_ret,
    int num_calls)
{
    (void)num_calls;
    TEST_ASSERT_EQUAL_PTR(make_context(0), context);
    TEST_ASSERT_EQUAL_PTR(make_device_id(0), device);
    TEST_ASSERT(properties[0] == CL_QUEUE_PROPERTIES);
    TEST_ASSERT(properties[1] == 0);
    return make_command_queue(0);
}
#endif // #if CL_HPP_TARGET_OPENCL_VERSION >= 200

void testCommandQueueFromSpecifiedContext()
{
    cl_command_queue expected = make_command_queue(0);
    cl_context expected_context =  make_context(0);
    cl_device_id expected_device = make_device_id(0);

    int context_refcount = 1;
    int device_refcount = 1;
    prepare_contextRefcounts(1, &expected_context, &context_refcount);
    prepare_deviceRefcounts(1, &expected_device, &device_refcount);

    // This is the context we will pass in to test
    cl::Context context = contextPool[0];

    // Assumes the context contains the fi rst device
    clGetContextInfo_StubWithCallback(clGetContextInfo_device);

#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clCreateCommandQueueWithProperties_StubWithCallback(clCreateCommandQueueWithProperties_testCommandQueueFromSpecifiedContext);
#else // #if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clCreateCommandQueue_StubWithCallback(clCreateCommandQueue_testCommandQueueFromSpecifiedContext);
#endif // #if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);

#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
#else // #if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
#endif // #if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clReleaseCommandQueue_ExpectAndReturn(expected, CL_SUCCESS);

    cl::CommandQueue queue(context);
    TEST_ASSERT_EQUAL_PTR(expected, queue());

    // Context not destroyed yet
    TEST_ASSERT_EQUAL(2, context_refcount);
    // Device object destroyed at end of scope
    TEST_ASSERT_EQUAL(1, device_refcount);

}

/****************************************************************************
 * Tests for cl::Device
 ****************************************************************************/

void testCopyDeviceNonNull1_1()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);

    cl::Device d0(make_device_id(0));
    cl::Device d1(make_device_id(1));
    d0 = d1;
}

void testCopyDeviceNonNull1_2()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);
    clRetainDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);

    cl::Device d0(make_device_id(0));
    cl::Device d1(make_device_id(1));
    d0 = d1;

    // Prevent destructor from interfering with the test
    d0() = NULL;
    d1() = NULL;
}

void testCopyDeviceFromNull1_1()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);
    // No other calls expected

    cl::Device d(make_device_id(0));
    d = cl::Device();
}

void testCopyDeviceFromNull1_2()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    cl::Device d(make_device_id(0));
    d = cl::Device();
}

void testCopyDeviceToNull1_1()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);
    // No other calls expected

    cl::Device d0;
    cl::Device d1(make_device_id(0));
    d0 = d1;
}

void testCopyDeviceToNull1_2()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clRetainDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    cl::Device d0;
    cl::Device d1(make_device_id(0));
    d0 = d1;

    // Prevent destructor from interfering with the test
    d0() = NULL;
    d1() = NULL;
}

void testCopyDeviceSelf()
{
    // Use 1.2 to check the retain/release calls
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);
    clRetainDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);

    cl::Device d0(make_device_id(0));
    cl::Device d1(make_device_id(1));
    d0 = d1;

    // Prevent destructor from interfering with the test
    d0() = NULL;
    d1() = NULL;
}

void testAssignDeviceNull()
{
    // Any version will do here
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    cl::Device d(make_device_id(0));
    d = (cl_device_id) NULL;
}

// These tests do not use the MAKE_MOVE_TESTS helper because they need to
// check whether the device is reference-countable, and to check that
// the reference-countable flag is correctly moved.
void testMoveAssignDeviceNonNull()
{
#ifdef TEST_RVALUE_REFERENCES
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);

    // Release called when trg overwritten
    clReleaseDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);

    cl::Device src(make_device_id(0));
    cl::Device trg(make_device_id(1));
    trg = std::move(src);
    TEST_ASSERT_EQUAL_PTR(make_device_id(0), trg());
    TEST_ASSERT_NULL(src());

    // Prevent destructor from interfering with the test
    trg() = NULL;
#endif
}

void testMoveAssignDeviceNull()
{
#ifdef TEST_RVALUE_REFERENCES
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);

    cl::Device trg;
    cl::Device src(make_device_id(1));
    trg = std::move(src);
    TEST_ASSERT_EQUAL_PTR(make_device_id(1), trg());
    TEST_ASSERT_NULL(src());

    // Prevent destructor from interfering with the test
    trg() = NULL;
#endif
}

void testMoveConstructDeviceNonNull()
{
#ifdef TEST_RVALUE_REFERENCES
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);

    cl::Device src(make_device_id(0));
    cl::Device trg(std::move(src));
    TEST_ASSERT_EQUAL_PTR(make_device_id(0), trg());
    TEST_ASSERT_NULL(src());

    // Prevent destructor from interfering with the test
    trg() = NULL;
#endif
}

void testMoveConstructDeviceNull()
{
#ifdef TEST_RVALUE_REFERENCES
    cl::Device empty;
    cl::Device trg(std::move(empty));
    TEST_ASSERT_NULL(trg());
    TEST_ASSERT_NULL(empty());
#endif
}

void testDestroyDevice1_1()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);
    // No other calls expected

    cl::Device d(make_device_id(0));
}

void testDestroyDevice1_2()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    cl::Device d(make_device_id(0));
}

static cl_int clGetDeviceIDs_PlatformWithZeroDevices(
    cl_platform_id  platform,
    cl_device_type  device_type,
    cl_uint  num_entries,
    cl_device_id  *devices,
    cl_uint  *num_devices,
    int num_calls)
{
    if (num_calls == 0)
    {
        TEST_ASSERT_EQUAL_PTR(make_platform_id(0), platform);
        TEST_ASSERT_EQUAL(CL_DEVICE_TYPE_ALL, device_type);
        TEST_ASSERT_NOT_NULL(num_devices);
        return CL_DEVICE_NOT_FOUND;
    }
    else
    {
        TEST_FAIL_MESSAGE("clGetDeviceIDs called too many times");
        return CL_INVALID_VALUE;
    }
}

void testPlatformWithZeroDevices()
{
    clGetDeviceIDs_StubWithCallback(clGetDeviceIDs_PlatformWithZeroDevices);

    cl::Platform p(make_platform_id(0));
    std::vector<cl::Device> devices;

    cl_int errCode = p.getDevices(CL_DEVICE_TYPE_ALL, &devices);
    TEST_ASSERT_EQUAL(CL_SUCCESS, errCode);
    TEST_ASSERT_EQUAL(0, devices.size());
}

/****************************************************************************
 * Tests for cl::Buffer
 ****************************************************************************/

void testMoveAssignBufferNonNull();
void testMoveAssignBufferNull();
void testMoveConstructBufferNonNull();
void testMoveConstructBufferNull();
MAKE_MOVE_TESTS(Buffer, make_mem, clReleaseMemObject, bufferPool);

// Stub of clCreateBuffer for testBufferConstructorContextInterator
// - return the first memory location

static cl_mem clCreateBuffer_testBufferConstructorContextIterator(
    cl_context context,
    cl_mem_flags flags,
    size_t size,
    void *host_ptr,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL_PTR(make_context(0), context);
    TEST_ASSERT_BITS(CL_MEM_COPY_HOST_PTR, flags, !CL_MEM_COPY_HOST_PTR);
    TEST_ASSERT_BITS(CL_MEM_READ_ONLY, flags, CL_MEM_READ_ONLY);
    TEST_ASSERT_EQUAL(sizeof(int)*1024, size);
    TEST_ASSERT_NULL(host_ptr);
    if (errcode_ret)
        errcode_ret = CL_SUCCESS;
    return make_mem(0);
}

// Declare forward these functions
static void * clEnqueueMapBuffer_testCopyHostToBuffer(
    cl_command_queue command_queue,
    cl_mem buffer,
    cl_bool blocking_map,
    cl_map_flags map_flags,
    size_t offset,
    size_t size,
    cl_uint num_events_in_wait_list,
    const cl_event *event_wait_list,
    cl_event *event,
    cl_int *errcode_ret,
    int num_calls);

static cl_int clEnqueueUnmapMemObject_testCopyHostToBuffer(
    cl_command_queue  command_queue ,
    cl_mem  memobj,
    void  *mapped_ptr,
    cl_uint  num_events_in_wait_list ,
    const cl_event  *event_wait_list ,
    cl_event  *event,
    int num_calls);

static cl_int clWaitForEvents_testCopyHostToBuffer(
    cl_uint num_events,
    const cl_event *event_list,
    int num_calls);

static cl_int clReleaseEvent_testCopyHostToBuffer(
    cl_event event,
    int num_calls);

void testBufferConstructorContextIterator()
{
    cl_mem expected = make_mem(0);

    // Assume this context includes make_device_id(0) for stub clGetContextInfo_device
    cl::Context context(make_context(0));

    clCreateBuffer_StubWithCallback(clCreateBuffer_testBufferConstructorContextIterator);
    clGetContextInfo_StubWithCallback(clGetContextInfo_device);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
#else // #if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
#endif //#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clRetainDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clCreateCommandQueueWithProperties_StubWithCallback(clCreateCommandQueueWithProperties_testCommandQueueFromSpecifiedContext);
#else // #if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clCreateCommandQueue_StubWithCallback(clCreateCommandQueue_testCommandQueueFromSpecifiedContext);
#endif //#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);
    clEnqueueMapBuffer_StubWithCallback(clEnqueueMapBuffer_testCopyHostToBuffer);
    clEnqueueUnmapMemObject_StubWithCallback(clEnqueueUnmapMemObject_testCopyHostToBuffer);
    clWaitForEvents_StubWithCallback(clWaitForEvents_testCopyHostToBuffer);
    clReleaseEvent_StubWithCallback(clReleaseEvent_testCopyHostToBuffer);
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(0), CL_SUCCESS);

    std::vector<int> host(1024);

    cl::Buffer buffer(context, host.begin(), host.end(), true);

    TEST_ASSERT_EQUAL_PTR(expected, buffer());

    // Tidy up at end of test
    clReleaseMemObject_ExpectAndReturn(expected, CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);
}

void testBufferConstructorQueueIterator()
{
    cl_context expected_context = make_context(0);
    int context_refcount = 1;
    cl_mem expected = make_mem(0);

    cl::CommandQueue queue(make_command_queue(0));

    prepare_contextRefcounts(1, &expected_context, &context_refcount);
    clGetCommandQueueInfo_StubWithCallback(clGetCommandQueueInfo_context);
    clCreateBuffer_StubWithCallback(clCreateBuffer_testBufferConstructorContextIterator);

    clEnqueueMapBuffer_StubWithCallback(clEnqueueMapBuffer_testCopyHostToBuffer);
    clEnqueueUnmapMemObject_StubWithCallback(clEnqueueUnmapMemObject_testCopyHostToBuffer);
    clWaitForEvents_StubWithCallback(clWaitForEvents_testCopyHostToBuffer);
    clReleaseEvent_StubWithCallback(clReleaseEvent_testCopyHostToBuffer);

    std::vector<int> host(1024);

    cl::Buffer buffer(queue, host.begin(), host.end(), true);

    TEST_ASSERT_EQUAL_PTR(expected, buffer());
    TEST_ASSERT_EQUAL(1, context_refcount);

    // Tidy up at end of test
    clReleaseMemObject_ExpectAndReturn(expected, CL_SUCCESS);
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(0), CL_SUCCESS);
}

/****************************************************************************
 * Tests for cl::Image1DBuffer
 ****************************************************************************/

/**
 * Stub for querying CL_IMAGE_BUFFER and returning make_mem(1).
 */
cl_int clGetImageInfo_testGetImageInfoBuffer(
    cl_mem image, cl_image_info param_name,
    size_t param_value_size, void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(0, num_calls);
    TEST_ASSERT_EQUAL_PTR(make_mem(0), image);
    TEST_ASSERT_EQUAL_HEX(CL_IMAGE_BUFFER, param_name);
    TEST_ASSERT_EQUAL(sizeof(cl_mem), param_value_size);

    if (param_value != NULL)
    {
        *(cl_mem *) param_value = make_mem(1);
    }
    if (param_value_size_ret != NULL)
        *param_value_size_ret = sizeof(cl_mem);
    return CL_SUCCESS;
}

void testGetImageInfoBuffer()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    cl_mem expected = make_mem(1);
    int refcount = 1;

    clGetImageInfo_StubWithCallback(clGetImageInfo_testGetImageInfoBuffer);
    prepare_memRefcounts(1, &expected, &refcount);

    cl::Image1DBuffer image(make_mem(0));
    const cl::Buffer &buffer = image.getImageInfo<CL_IMAGE_BUFFER>();
    TEST_ASSERT_EQUAL_PTR(make_mem(1), buffer());
    // Ref count should be 2 here because buffer has not been destroyed yet
    TEST_ASSERT_EQUAL(2, refcount);

    // prevent destructor from interfering with the test
    image() = NULL;
#endif
}

/**
 * Stub for querying CL_IMAGE_BUFFER and returning NULL.
 */
cl_int clGetImageInfo_testGetImageInfoBufferNull(
    cl_mem image, cl_image_info param_name,
    size_t param_value_size, void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(0, num_calls);
    TEST_ASSERT_EQUAL_PTR(make_mem(0), image);
    TEST_ASSERT_EQUAL_HEX(CL_IMAGE_BUFFER, param_name);
    TEST_ASSERT_EQUAL(sizeof(cl_mem), param_value_size);

    if (param_value != NULL)
    {
        *(cl_mem *) param_value = NULL;
    }
    if (param_value_size_ret != NULL)
        *param_value_size_ret = sizeof(cl_mem);
    return CL_SUCCESS;
}

void testGetImageInfoBufferNull()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clGetImageInfo_StubWithCallback(clGetImageInfo_testGetImageInfoBufferNull);

    cl::Image2D image(make_mem(0));
    cl::Buffer buffer = image.getImageInfo<CL_IMAGE_BUFFER>();
    TEST_ASSERT_NULL(buffer());

    // prevent destructor from interfering with the test
    image() = NULL;
#endif
}

void testGetImageInfoBufferOverwrite()
{
    clGetImageInfo_StubWithCallback(clGetImageInfo_testGetImageInfoBuffer);
    clReleaseMemObject_ExpectAndReturn(make_mem(2), CL_SUCCESS);
    clRetainMemObject_ExpectAndReturn(make_mem(1), CL_SUCCESS);

    cl::Image2D image(make_mem(0));
    cl::Buffer buffer(make_mem(2));
    cl_int status = image.getImageInfo(CL_IMAGE_BUFFER, &buffer);
    TEST_ASSERT_EQUAL(CL_SUCCESS, status);
    TEST_ASSERT_EQUAL_PTR(make_mem(1), buffer());

    // prevent destructor from interfering with the test
    image() = NULL;
    buffer() = NULL;
}

/**
 * A stub for clCreateImage that creates an image from a buffer
 * passing the buffer's cl_mem straight through.
 */
cl_mem clCreateImage_image1dbuffer(
    cl_context context,
    cl_mem_flags flags,
    const cl_image_format *image_format,
    const cl_image_desc *image_desc,
    void *host_ptr,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_NOT_NULL(image_format);
    TEST_ASSERT_NOT_NULL(image_desc);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_OBJECT_IMAGE1D_BUFFER, image_desc->image_type);

    // Return the passed buffer as the cl_mem
    return image_desc->buffer;
}

void testConstructImageFromBuffer()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 120
    const size_t width = 64;
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clCreateImage_StubWithCallback(clCreateImage_image1dbuffer);
    clReleaseMemObject_ExpectAndReturn(make_mem(0), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);

    cl::Context context(make_context(0));
    cl::Buffer buffer(make_mem(0));
    cl::Image1DBuffer image(
        context,
        CL_MEM_READ_ONLY,
        cl::ImageFormat(CL_R, CL_SIGNED_INT32),
        width,
        buffer);

    // Check that returned buffer matches the original
    TEST_ASSERT_EQUAL_PTR(buffer(), image());

    buffer() = NULL;
#endif
}

/****************************************************************************
 * Tests for cl::Image2D
 ****************************************************************************/

void testMoveAssignImage2DNonNull();
void testMoveAssignImage2DNull();
void testMoveConstructImage2DNonNull();
void testMoveConstructImage2DNull();
MAKE_MOVE_TESTS(Image2D, make_mem, clReleaseMemObject, image2DPool);

#ifdef CL_USE_DEPRECATED_OPENCL_1_1_APIS
static cl_mem clCreateImage2D_testCreateImage2D_1_1(
    cl_context context,
    cl_mem_flags flags,
    const cl_image_format *image_format,
    size_t image_width,
    size_t image_height,
    size_t image_row_pitch,
    void *host_ptr,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(0, num_calls);
    TEST_ASSERT_EQUAL_PTR(make_context(0), context);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_READ_WRITE, flags);

    TEST_ASSERT_NOT_NULL(image_format);
    TEST_ASSERT_EQUAL_HEX(CL_R, image_format->image_channel_order);
    TEST_ASSERT_EQUAL_HEX(CL_FLOAT, image_format->image_channel_data_type);

    TEST_ASSERT_EQUAL(64, image_width);
    TEST_ASSERT_EQUAL(32, image_height);
    TEST_ASSERT_EQUAL(256, image_row_pitch);
    TEST_ASSERT_NULL(host_ptr);

    if (errcode_ret != NULL)
        *errcode_ret = CL_SUCCESS;
    return make_mem(0);
}
#endif

void testCreateImage2D_1_1()
{
#ifdef CL_USE_DEPRECATED_OPENCL_1_1_APIS
    clGetContextInfo_StubWithCallback(clGetContextInfo_device);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);
    clCreateImage2D_StubWithCallback(clCreateImage2D_testCreateImage2D_1_1);

    cl_int err;
    cl::Context context;
    context() = make_context(0);
    cl::Image2D image(
        context, CL_MEM_READ_WRITE,
        cl::ImageFormat(CL_R, CL_FLOAT), 64, 32, 256, NULL, &err);

    TEST_ASSERT_EQUAL(CL_SUCCESS, err);
    TEST_ASSERT_EQUAL_PTR(make_mem(0), image());

    context() = NULL;
    image() = NULL;
#endif
}

static cl_mem clCreateImage_testCreateImage2D_1_2(
    cl_context context,
    cl_mem_flags flags,
    const cl_image_format *image_format,
    const cl_image_desc *image_desc,
    void *host_ptr,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(0, num_calls);
    TEST_ASSERT_EQUAL_PTR(make_context(0), context);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_READ_WRITE, flags);

    TEST_ASSERT_NOT_NULL(image_format);
    TEST_ASSERT_EQUAL_HEX(CL_R, image_format->image_channel_order);
    TEST_ASSERT_EQUAL_HEX(CL_FLOAT, image_format->image_channel_data_type);

    TEST_ASSERT_NOT_NULL(image_desc);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_OBJECT_IMAGE2D, image_desc->image_type);
    TEST_ASSERT_EQUAL(64, image_desc->image_width);
    TEST_ASSERT_EQUAL(32, image_desc->image_height);
    TEST_ASSERT_EQUAL(256, image_desc->image_row_pitch);
    TEST_ASSERT_EQUAL(0, image_desc->num_mip_levels);
    TEST_ASSERT_EQUAL(0, image_desc->num_samples);
    TEST_ASSERT_NULL(image_desc->buffer);

    TEST_ASSERT_NULL(host_ptr);

    if (errcode_ret != NULL)
        *errcode_ret = CL_SUCCESS;
    return make_mem(0);
}

void testCreateImage2D_1_2()
{
    clGetContextInfo_StubWithCallback(clGetContextInfo_device);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clCreateImage_StubWithCallback(clCreateImage_testCreateImage2D_1_2);

    cl_int err;
    cl::Context context;
    context() = make_context(0);
    cl::Image2D image(
        context, CL_MEM_READ_WRITE,
        cl::ImageFormat(CL_R, CL_FLOAT), 64, 32, 256, NULL, &err);

    TEST_ASSERT_EQUAL(CL_SUCCESS, err);
    TEST_ASSERT_EQUAL_PTR(make_mem(0), image());

    context() = NULL;
    image() = NULL;
}

/****************************************************************************
 * Tests for cl::Image3D
 ****************************************************************************/

void testMoveAssignImage3DNonNull();
void testMoveAssignImage3DNull();
void testMoveConstructImage3DNonNull();
void testMoveConstructImage3DNull();
MAKE_MOVE_TESTS(Image3D, make_mem, clReleaseMemObject, image3DPool);

#ifdef CL_USE_DEPRECATED_OPENCL_1_1_APIS
static cl_mem clCreateImage3D_testCreateImage3D_1_1(
    cl_context context,
    cl_mem_flags flags,
    const cl_image_format *image_format,
    size_t image_width,
    size_t image_height,
    size_t image_depth,
    size_t image_row_pitch,
    size_t image_slice_pitch,
    void *host_ptr,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(0, num_calls);
    TEST_ASSERT_EQUAL_PTR(make_context(0), context);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, flags);

    TEST_ASSERT_NOT_NULL(image_format);
    TEST_ASSERT_EQUAL_HEX(CL_R, image_format->image_channel_order);
    TEST_ASSERT_EQUAL_HEX(CL_FLOAT, image_format->image_channel_data_type);

    TEST_ASSERT_EQUAL(64, image_width);
    TEST_ASSERT_EQUAL(32, image_height);
    TEST_ASSERT_EQUAL(16, image_depth);
    TEST_ASSERT_EQUAL(256, image_row_pitch);
    TEST_ASSERT_EQUAL(65536, image_slice_pitch);
    TEST_ASSERT_EQUAL_PTR((void *) 0xdeadbeef, host_ptr);

    if (errcode_ret != NULL)
        *errcode_ret = CL_SUCCESS;
    return make_mem(0);
}
#endif

void testCreateImage3D_1_1()
{
#ifdef CL_USE_DEPRECATED_OPENCL_1_1_APIS
    clGetContextInfo_StubWithCallback(clGetContextInfo_device);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_1);
    clCreateImage3D_StubWithCallback(clCreateImage3D_testCreateImage3D_1_1);

    cl_int err;
    cl::Context context;
    context() = make_context(0);
    cl::Image3D image(
        context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        cl::ImageFormat(CL_R, CL_FLOAT), 64, 32, 16, 256, 65536, (void *) 0xdeadbeef, &err);

    TEST_ASSERT_EQUAL(CL_SUCCESS, err);
    TEST_ASSERT_EQUAL_PTR(make_mem(0), image());

    context() = NULL;
    image() = NULL;
#endif
}

static cl_mem clCreateImage_testCreateImage3D_1_2(
    cl_context context,
    cl_mem_flags flags,
    const cl_image_format *image_format,
    const cl_image_desc *image_desc,
    void *host_ptr,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(0, num_calls);
    TEST_ASSERT_EQUAL_PTR(make_context(0), context);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, flags);

    TEST_ASSERT_NOT_NULL(image_format);
    TEST_ASSERT_EQUAL_HEX(CL_R, image_format->image_channel_order);
    TEST_ASSERT_EQUAL_HEX(CL_FLOAT, image_format->image_channel_data_type);

    TEST_ASSERT_NOT_NULL(image_desc);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_OBJECT_IMAGE3D, image_desc->image_type);
    TEST_ASSERT_EQUAL(64, image_desc->image_width);
    TEST_ASSERT_EQUAL(32, image_desc->image_height);
    TEST_ASSERT_EQUAL(16, image_desc->image_depth);
    TEST_ASSERT_EQUAL(256, image_desc->image_row_pitch);
    TEST_ASSERT_EQUAL(65536, image_desc->image_slice_pitch);
    TEST_ASSERT_EQUAL(0, image_desc->num_mip_levels);
    TEST_ASSERT_EQUAL(0, image_desc->num_samples);
    TEST_ASSERT_NULL(image_desc->buffer);

    TEST_ASSERT_EQUAL_PTR((void *) 0xdeadbeef, host_ptr);

    if (errcode_ret != NULL)
        *errcode_ret = CL_SUCCESS;
    return make_mem(0);
}

void testCreateImage3D_1_2()
{
    clGetContextInfo_StubWithCallback(clGetContextInfo_device);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clCreateImage_StubWithCallback(clCreateImage_testCreateImage3D_1_2);

    cl_int err;
    cl::Context context;
    context() = make_context(0);
    cl::Image3D image(
        context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        cl::ImageFormat(CL_R, CL_FLOAT), 64, 32, 16, 256, 65536, (void *) 0xdeadbeef, &err);

    TEST_ASSERT_EQUAL(CL_SUCCESS, err);
    TEST_ASSERT_EQUAL_PTR(make_mem(0), image());

    context() = NULL;
    image() = NULL;
}

/****************************************************************************
 * Tests for cl::Kernel
 ****************************************************************************/
void testMoveAssignKernelNonNull();
void testMoveAssignKernelNull();
void testMoveConstructKernelNonNull();
void testMoveConstructKernelNull();
MAKE_MOVE_TESTS(Kernel, make_kernel, clReleaseKernel, kernelPool);

static cl_int scalarArg;
static cl_int3 vectorArg;

void testKernelSetArgScalar()
{
    scalarArg = 0xcafebabe;
    clSetKernelArg_ExpectAndReturn(make_kernel(0), 3, 4, &scalarArg, CL_SUCCESS);
    kernelPool[0].setArg(3, scalarArg);
}

void testKernelSetArgVector()
{
    vectorArg.s[0] = 0x12345678;
    vectorArg.s[1] = 0x23456789;
    vectorArg.s[2] = 0x87654321;
    clSetKernelArg_ExpectAndReturn(make_kernel(0), 2, 16, &vectorArg, CL_SUCCESS);
    kernelPool[0].setArg(2, vectorArg);
}

void testKernelSetArgMem()
{
    clSetKernelArg_ExpectAndReturn(make_kernel(0), 1, sizeof(cl_mem), &bufferPool[1](), CL_SUCCESS);
    kernelPool[0].setArg(1, bufferPool[1]);
}

void testKernelSetArgLocal()
{
    clSetKernelArg_ExpectAndReturn(make_kernel(0), 2, 123, NULL, CL_SUCCESS);
    kernelPool[0].setArg(2, cl::Local(123));
}

void testKernelSetExecInfo()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    cl_bool val = CL_TRUE;
    // Using CL_KERNEL_EXEC_INFO_SVM_FINE_GRAIN_SYSTEM in the tests since it's
    // defined by the core spec but this function is particularly useful for
    // vendor extensions.
    clSetKernelExecInfo_ExpectAndReturn(make_kernel(0),
                                        CL_KERNEL_EXEC_INFO_SVM_FINE_GRAIN_SYSTEM,
                                        sizeof(cl_bool), &val, CL_SUCCESS);
    kernelPool[0].setExecInfo(CL_KERNEL_EXEC_INFO_SVM_FINE_GRAIN_SYSTEM, val);
    // Also test the typesafe version
    clSetKernelExecInfo_ExpectAndReturn(make_kernel(0),
                                        CL_KERNEL_EXEC_INFO_SVM_FINE_GRAIN_SYSTEM,
                                        sizeof(cl_bool), &val, CL_SUCCESS);
    kernelPool[0].setExecInfo<CL_KERNEL_EXEC_INFO_SVM_FINE_GRAIN_SYSTEM>(val);
#endif
}

/****************************************************************************
 * Tests for cl::copy
 ****************************************************************************/

// This method should allocate some host accesible memory
// so we must do this ourselves
void *some_host_memory;

static void * clEnqueueMapBuffer_testCopyHostToBuffer(
    cl_command_queue command_queue,
    cl_mem buffer,
    cl_bool blocking_map,
    cl_map_flags map_flags,
    size_t offset,
    size_t size,
    cl_uint num_events_in_wait_list,
    const cl_event *event_wait_list,
    cl_event *event,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL_PTR(make_command_queue(0), command_queue);
    TEST_ASSERT_EQUAL_PTR(make_mem(0), buffer);
    TEST_ASSERT_EQUAL(CL_TRUE, blocking_map);
    TEST_ASSERT_EQUAL(CL_MAP_WRITE, map_flags);
    TEST_ASSERT_EQUAL(sizeof(int)*1024, size);

    some_host_memory = malloc(sizeof(int) * 1024);

    // Set the return event
    if (event)
        *event = NULL;

    // Set the return error code
    if (errcode_ret)
        *errcode_ret = CL_SUCCESS;

    return some_host_memory;
}

static cl_int clEnqueueUnmapMemObject_testCopyHostToBuffer(
    cl_command_queue  command_queue ,
    cl_mem  memobj,
    void  *mapped_ptr,
    cl_uint  num_events_in_wait_list ,
    const cl_event  *event_wait_list ,
    cl_event  *event,
    int num_calls)
{
    TEST_ASSERT_EQUAL_PTR(make_command_queue(0), command_queue);
    TEST_ASSERT_EQUAL_PTR(make_mem(0), memobj);
    TEST_ASSERT_EQUAL_PTR(some_host_memory, mapped_ptr);
    TEST_ASSERT_NOT_NULL(event);
    return CL_SUCCESS;
}

static cl_int clWaitForEvents_testCopyHostToBuffer(
    cl_uint num_events,
    const cl_event *event_list,
    int num_calls)
{
    TEST_ASSERT_NOT_NULL(event_list);
    TEST_ASSERT_EQUAL(1, num_events);
    return CL_SUCCESS;
}

static cl_int clReleaseEvent_testCopyHostToBuffer(
    cl_event event,
    int num_calls)
{
    TEST_ASSERT_NOT_NULL(event);
    return CL_SUCCESS;
}

void testCopyHostToBuffer()
{
    cl_context context_expect = make_context(0);
    int context_refcount = 1;
    prepare_contextRefcounts(1, &context_expect, &context_refcount);
    cl::Context context = contextPool[0];

    cl_mem mem_expect = make_mem(0);
    int mem_refcount = 1;
    prepare_memRefcounts(1, &mem_expect, &mem_refcount);
    cl::Buffer buffer(make_mem(0));

    cl_command_queue queue_expect = make_command_queue(0);
    cl::CommandQueue queue(queue_expect);
    clReleaseCommandQueue_ExpectAndReturn(queue_expect, CL_SUCCESS);

    // Returns the pointer to host memory
    clEnqueueMapBuffer_StubWithCallback(clEnqueueMapBuffer_testCopyHostToBuffer);
    clEnqueueUnmapMemObject_StubWithCallback(clEnqueueUnmapMemObject_testCopyHostToBuffer);

    clWaitForEvents_StubWithCallback(clWaitForEvents_testCopyHostToBuffer);
    clReleaseEvent_StubWithCallback(clReleaseEvent_testCopyHostToBuffer);

    std::vector<int> host(1024);
    for (int i = 0; i < 1024; i++)
        host[i] = i;

    cl::copy(queue, host.begin(), host.end(), buffer);

    // Check that the memory was copied to some_host_memory
    TEST_ASSERT_EQUAL_MEMORY(&host[0], some_host_memory, sizeof(int) * 1024);

    free(some_host_memory);

}

/****************************************************************************
* Tests for building Programs
****************************************************************************/

static cl_int clGetDeviceInfo_testGetBuildInfo(
    cl_device_id device,
    cl_device_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(param_name, CL_DEVICE_PLATFORM);
    TEST_ASSERT_EQUAL(param_value_size, sizeof(cl_platform_id));
    TEST_ASSERT_NOT_EQUAL(param_value, NULL);
    TEST_ASSERT_EQUAL(param_value_size_ret, NULL);
    cl_platform_id temp = make_platform_id(0);
    memcpy(param_value, &temp, sizeof(cl_platform_id));
    return CL_SUCCESS;
}


static  cl_int clGetProgramBuildInfo_testGetBuildInfo(
    cl_program program,
    cl_device_id device,
    cl_program_build_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(param_name, CL_PROGRAM_BUILD_LOG);

    const char returnString[] = 
        "This is the string returned by the build info function.";
    if (param_value) {
        ::size_t returnSize = param_value_size;
        if (sizeof(returnString) < returnSize) {
            returnSize = sizeof(returnString);
        }
        memcpy(param_value, returnString, returnSize);
    }
    else {
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(returnString);
        }
    }

    return CL_SUCCESS;
}

void testGetBuildInfo()
{
    cl_device_id fakeDevice = make_device_id(0);
    clGetDeviceInfo_ExpectAndReturn(fakeDevice, CL_DEVICE_PLATFORM, sizeof(cl_platform_id), NULL, NULL, CL_SUCCESS);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_testGetBuildInfo);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);
    clGetProgramBuildInfo_StubWithCallback(clGetProgramBuildInfo_testGetBuildInfo);
    clGetProgramBuildInfo_StubWithCallback(clGetProgramBuildInfo_testGetBuildInfo);

    cl::Program prog(make_program(0));
    cl::Device dev(fakeDevice);
    
    cl_int err;
    std::string log = prog.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev, &err);

    prog() = NULL;
    dev() = NULL;
}

static cl_int clBuildProgram_testBuildProgram(
    cl_program           program,
    cl_uint              num_devices,
    const cl_device_id * device_list,
    const char *         options,
    void (CL_CALLBACK *  pfn_notify)(cl_program program, void * user_data),
    void *               user_data,
    int num_calls)
{
    TEST_ASSERT_EQUAL(program, make_program(0));
    TEST_ASSERT_NOT_EQUAL(num_devices, 0);
    TEST_ASSERT_NOT_EQUAL(device_list, NULL);
    TEST_ASSERT_EQUAL(options, NULL);
    TEST_ASSERT_EQUAL(pfn_notify, NULL);
    TEST_ASSERT_EQUAL(user_data, NULL);

    for (cl_uint i = 0; i < num_devices; i++) {
        TEST_ASSERT_EQUAL(device_list[i], make_device_id(i));
    }

    return CL_SUCCESS;
}

void testBuildProgramSingleDevice()
{
    cl_program program = make_program(0);
    cl_device_id device_id = make_device_id(0);
    int sc = 0;

    // Creating a device queries the platform version:
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_1_2);

    clBuildProgram_StubWithCallback(clBuildProgram_testBuildProgram);

    // Building the program queries the program build log:
    clRetainDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);
    clGetProgramBuildInfo_StubWithCallback(clGetProgramBuildInfo_testGetBuildInfo);
    clGetProgramBuildInfo_StubWithCallback(clGetProgramBuildInfo_testGetBuildInfo);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    clReleaseProgram_ExpectAndReturn(program, CL_SUCCESS);

    cl::Program prog(program);
    cl::Device dev(device_id);

    cl_int errcode = prog.build(dev);

    TEST_ASSERT_EQUAL(errcode, CL_SUCCESS);
}

/**
* Stub implementation of clGetCommandQueueInfo that returns first one image then none
*/
static cl_int clGetSupportedImageFormats_testGetSupportedImageFormats(
    cl_context context,
    cl_mem_flags flags,
    cl_mem_object_type image_type,
    cl_uint num_entries,
    cl_image_format *image_formats,
    cl_uint *num_image_formats,
    int num_calls)
{        
    // Catch failure case that causes error in bugzilla 13355:
    // returns CL_INVALID_VALUE if flags or image_type are not valid, 
    // or if num_entries is 0 and image_formats is not NULL.
    if (num_entries == 0 && image_formats != NULL) {
        return CL_INVALID_VALUE;
    }
    if (num_entries == 0)  {
        // If num_entries was 0 this is the query for number
        if (num_image_formats) {
            if (num_calls == 0) {
                *num_image_formats = 1;
            }
            else {
                *num_image_formats = 0;
            }
        }
    }
    else {
        // Should return something
        TEST_ASSERT_NOT_NULL(image_formats);
        
        // For first call we should return one format here
        if (num_calls == 1) {
            TEST_ASSERT_EQUAL(num_entries, 1);
            image_formats[0] = cl::ImageFormat(CL_RGB, CL_FLOAT);
        }
    }

    return CL_SUCCESS;
}

void testGetSupportedImageFormats()
{
    cl_context ctx_cl = make_context(0);

    clGetSupportedImageFormats_StubWithCallback(clGetSupportedImageFormats_testGetSupportedImageFormats);
    clGetSupportedImageFormats_StubWithCallback(clGetSupportedImageFormats_testGetSupportedImageFormats);
    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);

    cl::Context ctx(ctx_cl);
    std::vector<cl::ImageFormat> formats;
    cl_int ret = CL_SUCCESS;

    ret = ctx.getSupportedImageFormats(
        CL_MEM_READ_WRITE,
        CL_MEM_OBJECT_IMAGE2D,
        &formats);
    TEST_ASSERT_EQUAL(ret, CL_SUCCESS);
    TEST_ASSERT_EQUAL(formats.size(), 1);
    ret = ctx.getSupportedImageFormats(
        CL_MEM_READ_WRITE,
        CL_MEM_OBJECT_IMAGE2D,
        &formats);
    TEST_ASSERT_EQUAL(formats.size(), 0);
    TEST_ASSERT_EQUAL(ret, CL_SUCCESS);
}

void testCreateSubDevice()
{
    // TODO

}

void testGetContextInfoDevices()
{
    // TODO
}

static cl_mem clCreateImage_testCreateImage2DFromBuffer_2_0(
    cl_context context,
    cl_mem_flags flags,
    const cl_image_format *image_format,
    const cl_image_desc *image_desc,
    void *host_ptr,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_NOT_NULL(image_format);
    TEST_ASSERT_NOT_NULL(image_desc);
    TEST_ASSERT_NULL(host_ptr);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_OBJECT_IMAGE2D, image_desc->image_type);

    // Return the passed buffer as the cl_mem and success for the error code
    if (errcode_ret) {
        *errcode_ret = CL_SUCCESS;
    }
    return image_desc->buffer;
}

void testCreateImage2DFromBuffer_2_0()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clGetContextInfo_StubWithCallback(clGetContextInfo_device);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
    clCreateImage_StubWithCallback(clCreateImage_testCreateImage2DFromBuffer_2_0);
    clReleaseMemObject_ExpectAndReturn(make_mem(0), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);

    cl_int err;
    cl::Context context(make_context(0));

    // Create buffer
    // Create image from buffer
    cl::Buffer buffer(make_mem(0));
    cl::Image2D imageFromBuffer(
        context,
        cl::ImageFormat(CL_R, CL_FLOAT), buffer, 64, 32, 256, &err);

    TEST_ASSERT_EQUAL_PTR(buffer(), imageFromBuffer());
    TEST_ASSERT_EQUAL(CL_SUCCESS, err);

    buffer() = NULL;
#endif
}

static cl_mem clCreateImage_testCreateImage2D_2_0(
    cl_context context,
    cl_mem_flags flags,
    const cl_image_format *image_format,
    const cl_image_desc *image_desc,
    void *host_ptr,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_EQUAL(0, num_calls);
    TEST_ASSERT_EQUAL_PTR(make_context(0), context);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_READ_WRITE, flags);

    TEST_ASSERT_NOT_NULL(image_format);
    TEST_ASSERT_EQUAL_HEX(CL_RGBA, image_format->image_channel_order);
    TEST_ASSERT_EQUAL_HEX(CL_FLOAT, image_format->image_channel_data_type);

    TEST_ASSERT_NOT_NULL(image_desc);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_OBJECT_IMAGE2D, image_desc->image_type);
    TEST_ASSERT_EQUAL(64, image_desc->image_width);
    TEST_ASSERT_EQUAL(32, image_desc->image_height);
    TEST_ASSERT_EQUAL(256, image_desc->image_row_pitch);
    TEST_ASSERT_EQUAL(0, image_desc->num_mip_levels);
    TEST_ASSERT_EQUAL(0, image_desc->num_samples);
    TEST_ASSERT_NULL(image_desc->buffer);

    TEST_ASSERT_NULL(host_ptr);

    if (errcode_ret != NULL)
        *errcode_ret = CL_SUCCESS;
    return make_mem(0);
}

static cl_mem clCreateImage_testCreateImage2DFromImage_2_0(
    cl_context context,
    cl_mem_flags flags,
    const cl_image_format *image_format,
    const cl_image_desc *image_desc,
    void *host_ptr,
    cl_int *errcode_ret,
    int num_calls)
{
    TEST_ASSERT_NOT_NULL(image_format);
    TEST_ASSERT_NOT_NULL(image_desc);
    TEST_ASSERT_NULL(host_ptr);
    TEST_ASSERT_EQUAL_HEX(CL_MEM_OBJECT_IMAGE2D, image_desc->image_type);

    // Return the passed buffer as the cl_mem and success for the error code
    if (errcode_ret) {
        *errcode_ret = CL_SUCCESS;
    }
    return image_desc->buffer;
}

static cl_int clGetImageInfo_testCreateImage2DFromImage_2_0(
    cl_mem image,
    cl_image_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    TEST_ASSERT_INT_WITHIN(6, 0, num_calls);
    return CL_SUCCESS;
}

void testCreateImage2DFromImage_2_0()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clGetContextInfo_StubWithCallback(clGetContextInfo_device);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
    clCreateImage_StubWithCallback(clCreateImage_testCreateImage2D_2_0);


    cl_int err;
    cl::Context context(make_context(0));

    // As in 1.2 2D image test, needed as source for image-from-image
    cl::Image2D image(
        context, CL_MEM_READ_WRITE,
        cl::ImageFormat(CL_RGBA, CL_FLOAT), 64, 32, 256, NULL, &err);

    TEST_ASSERT_EQUAL(CL_SUCCESS, err);
    TEST_ASSERT_EQUAL_PTR(make_mem(0), image());

    // Continue state for next phase
    clGetImageInfo_StubWithCallback(clGetImageInfo_testCreateImage2DFromImage_2_0);
    clCreateImage_StubWithCallback(clCreateImage_testCreateImage2DFromImage_2_0);
    clReleaseMemObject_ExpectAndReturn(make_mem(0), CL_SUCCESS);
    clReleaseMemObject_ExpectAndReturn(make_mem(0), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(0), CL_SUCCESS);

    // Create 2D image from 2D Image with a new channel order
    cl::Image2D imageFromImage(
        context,
        CL_sRGB,
        image,
        &err
        );

    TEST_ASSERT_EQUAL(CL_SUCCESS, err);
    TEST_ASSERT_EQUAL_PTR(image(), imageFromImage());

    //imageFromImage() = NULL;
    //image() = NULL;
    //context() = NULL;
#endif
}

// Note that default tests maintain state when run from the same
// unit process.
// One default setting test will maintain the defaults until the end.
void testSetDefaultPlatform()
{
    cl::Platform p(make_platform_id(1));
    cl::Platform p2 = cl::Platform::setDefault(p);
    cl::Platform p3 = cl::Platform::getDefault();
    TEST_ASSERT_EQUAL(p(), p2());
    TEST_ASSERT_EQUAL(p(), p3());
}

// Note that default tests maintain state when run from the same
// unit process.
// One default setting test will maintain the defaults until the end.
void testSetDefaultPlatformTwice()
{
    cl::Platform p(make_platform_id(2));
    cl::Platform p2 = cl::Platform::getDefault();
    cl::Platform p3 = cl::Platform::setDefault(p);
    // Set default should have failed
    TEST_ASSERT_EQUAL(p2(), p3());
    TEST_ASSERT_NOT_EQUAL(p(), p3());
}

// Note that default tests maintain state when run from the same
// unit process.
// One default setting test will maintain the defaults until the end.
void testSetDefaultContext()
{   

    clRetainContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clRetainContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clRetainContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(1), CL_SUCCESS);

    cl::Context c(make_context(1));
    cl::Context c2 = cl::Context::setDefault(c);
    cl::Context c3 = cl::Context::getDefault();
    TEST_ASSERT_EQUAL(c(), c2());
    TEST_ASSERT_EQUAL(c(), c3());
}

// Note that default tests maintain state when run from the same
// unit process.
// One default setting test will maintain the defaults until the end.
void testSetDefaultCommandQueue()
{
    clRetainCommandQueue_ExpectAndReturn(make_command_queue(1), CL_SUCCESS);
    clRetainCommandQueue_ExpectAndReturn(make_command_queue(1), CL_SUCCESS);
    clRetainCommandQueue_ExpectAndReturn(make_command_queue(1), CL_SUCCESS);
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(1), CL_SUCCESS);
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(1), CL_SUCCESS);
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(1), CL_SUCCESS);

    cl::CommandQueue c(make_command_queue(1));
    cl::CommandQueue c2 = cl::CommandQueue::setDefault(c);
    cl::CommandQueue c3 = cl::CommandQueue::getDefault();
    TEST_ASSERT_EQUAL(c(), c2());
    TEST_ASSERT_EQUAL(c(), c3());
}

// Note that default tests maintain state when run from the same
// unit process.
// One default setting test will maintain the defaults until the end.
void testSetDefaultDevice()
{
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);

    clRetainDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);
    clRetainDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);
    clRetainDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);
    clReleaseDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);
    clReleaseDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);
    clReleaseDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);

    cl::Device d(make_device_id(1));
    cl::Device  d2 = cl::Device::setDefault(d);
    cl::Device  d3 = cl::Device::getDefault();
    TEST_ASSERT_EQUAL(d(), d2());
    TEST_ASSERT_EQUAL(d(), d3());
}

static cl_command_queue clCreateCommandQueueWithProperties_testCommandQueueDevice(
    cl_context context,
    cl_device_id device,
    const cl_queue_properties *properties,
    cl_int *errcode_ret,
    int num_calls)
{
    (void)num_calls;
    TEST_ASSERT_EQUAL_PTR(make_context(1), context);
    TEST_ASSERT_EQUAL_PTR(make_device_id(1), device);
    TEST_ASSERT_EQUAL(properties[0], CL_QUEUE_PROPERTIES);
    static cl_command_queue default_ = 0;

    if ((properties[1] & CL_QUEUE_ON_DEVICE_DEFAULT) == 0) {
        TEST_ASSERT_EQUAL(properties[1], (CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_ON_DEVICE));
        if (properties[2] == CL_QUEUE_SIZE) {
            TEST_ASSERT_EQUAL(properties[3], 256);
            TEST_ASSERT_EQUAL(properties[4], 0);
            return make_command_queue(2);
        }
        else {
            TEST_ASSERT_EQUAL(properties[2], 0);
            return make_command_queue(3);
        }
    }
    else {
        TEST_ASSERT_EQUAL(properties[1], (CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_ON_DEVICE | CL_QUEUE_ON_DEVICE_DEFAULT));
        if (default_ == 0) {
            default_ = make_command_queue(4);
        }
        return default_;
    }
}

void testCreateDeviceCommandQueue()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clRetainContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
    clCreateCommandQueueWithProperties_StubWithCallback(clCreateCommandQueueWithProperties_testCommandQueueDevice);       
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(4), CL_SUCCESS);
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(4), CL_SUCCESS);
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(2), CL_SUCCESS);
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(3), CL_SUCCESS);
    clReleaseDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(1), CL_SUCCESS);

    cl::Context c(make_context(1));
    cl::Context c2 = cl::Context::setDefault(c);
    cl::Device d(make_device_id(1));

    cl::DeviceCommandQueue dq(c, d);
    cl::DeviceCommandQueue dq2(c, d, 256);    

    cl::DeviceCommandQueue dqd = cl::DeviceCommandQueue::makeDefault(c, d);
    cl::DeviceCommandQueue dqd2 = cl::DeviceCommandQueue::makeDefault(c, d);

    TEST_ASSERT_EQUAL(dqd(), dqd2());
#endif
}

static cl_mem clCreatePipe_testCreatePipe(
    cl_context context,
    cl_mem_flags flags,
    cl_uint packet_size,
    cl_uint num_packets,
    const cl_pipe_properties *props,
    cl_int *errcode_ret,
    int num_calls)
{
    if (flags == 0) {
        flags = CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS;
    }
    TEST_ASSERT_EQUAL(flags, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS);
    TEST_ASSERT_NULL(props);

    if (errcode_ret)
        *errcode_ret = CL_SUCCESS;
    return make_mem(0);
}

static cl_int clGetPipeInfo_testCreatePipe(
    cl_mem pipe,
    cl_pipe_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    TEST_ASSERT_NOT_NULL(param_value);
    if (param_name == CL_PIPE_PACKET_SIZE) {
        *static_cast<cl_uint*>(param_value) = 16;
        if (param_value_size_ret) {
            *param_value_size_ret = param_value_size;
        }
        return CL_SUCCESS;
    }
    else if (param_name == CL_PIPE_MAX_PACKETS) {
        *static_cast<cl_uint*>(param_value) = 32;
        if (param_value_size_ret) {
            *param_value_size_ret = param_value_size;
        }
        return CL_SUCCESS;
    }
    else {
        TEST_FAIL();
        return CL_INVALID_VALUE;
    }
}

void testCreatePipe()
{    
#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    clCreatePipe_StubWithCallback(clCreatePipe_testCreatePipe);
    clGetPipeInfo_StubWithCallback(clGetPipeInfo_testCreatePipe);
    clRetainContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clReleaseMemObject_ExpectAndReturn(make_mem(0), CL_SUCCESS);
    clReleaseMemObject_ExpectAndReturn(make_mem(0), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(1), CL_SUCCESS);

    cl::Context c(make_context(1));
    cl::Pipe p(c, 16, 32);
    cl::Pipe p2(16, 32);

    cl_uint size = p2.getInfo<CL_PIPE_PACKET_SIZE>();
    cl_uint packets;
    p2.getInfo(CL_PIPE_MAX_PACKETS, &packets);

    TEST_ASSERT_EQUAL(size, 16);
    TEST_ASSERT_EQUAL(packets, 32);
#endif
}

static cl_int clGetKernelSubGroupInfo_testSubGroups(cl_kernel kernel,
    cl_device_id device,
    cl_kernel_sub_group_info param_name,
    size_t input_value_size,
    const void *input_value,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{    
    TEST_ASSERT_NOT_NULL(input_value);
    TEST_ASSERT_NOT_NULL(param_value);

    if (param_name == CL_KERNEL_MAX_SUB_GROUP_SIZE_FOR_NDRANGE_KHR) {
        *static_cast<size_t*>(param_value) = 32;
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(size_t);
        }
        return CL_SUCCESS;
    }
    else if (param_name == CL_KERNEL_SUB_GROUP_COUNT_FOR_NDRANGE_KHR) {
        *static_cast<size_t*>(param_value) = 2;
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(size_t);
        }
        return CL_SUCCESS;
    }
    else {
        TEST_ABORT();
        return CL_INVALID_OPERATION;
    }
}

void testSubGroups()
{
// TODO support testing cl_khr_subgroups on 2.0
#if CL_HPP_TARGET_OPENCL_VERSION >= 210
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
    clGetKernelSubGroupInfo_StubWithCallback(clGetKernelSubGroupInfo_testSubGroups);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);
    clReleaseKernel_ExpectAndReturn(make_kernel(0), CL_SUCCESS);

    cl::Kernel k(make_kernel(0));
    cl::Device d(make_device_id(0));
    cl_int err;
    cl::NDRange ndrange(8, 8);
    size_t res1 = k.getSubGroupInfo<CL_KERNEL_MAX_SUB_GROUP_SIZE_FOR_NDRANGE_KHR>(
        d, ndrange, &err);
    size_t res2 = 0;
    err = k.getSubGroupInfo(
        d, CL_KERNEL_SUB_GROUP_COUNT_FOR_NDRANGE_KHR, ndrange, &res2);

    TEST_ASSERT_EQUAL(res1, 32);
    TEST_ASSERT_EQUAL(res2, 2);
#endif
}

/**
* Stub implementation of clGetDeviceInfo that returns an absense of builtin kernels
*/
static cl_int clGetDeviceInfo_builtin(
    cl_device_id id,
    cl_device_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    // Test to verify case where empty string is returned - so size is 0
    (void)num_calls;
    TEST_ASSERT_EQUAL_HEX(CL_DEVICE_BUILT_IN_KERNELS, param_name);
    if (param_value == NULL) {
        if (param_value_size_ret != NULL) {
            *param_value_size_ret = 0;
        }
    }
    return CL_SUCCESS;
}

void testBuiltInKernels()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 120
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    cl::Device d0(make_device_id(0));

    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_builtin);
    cl::string s = d0.getInfo<CL_DEVICE_BUILT_IN_KERNELS>();
#endif
}

/**
 * Stub implementation of clCloneKernel that returns a new kernel object
 */
static cl_kernel clCloneKernel_simplecopy(
    cl_kernel k,
    cl_int *err,
    int num_calls)
{
    // Test to verify case where empty string is returned - so size is 0
    (void)num_calls;
    return make_kernel(POOL_MAX);
    return CL_SUCCESS;
}

void testCloneKernel()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 210
    clCloneKernel_StubWithCallback(clCloneKernel_simplecopy);
    clReleaseKernel_ExpectAndReturn(make_kernel(POOL_MAX), CL_SUCCESS);
    cl::Kernel clone = kernelPool[0].clone();
    TEST_ASSERT_EQUAL(clone(), make_kernel(POOL_MAX));
#endif
}

void testEnqueueMapSVM()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 200
    std::vector<int> vec(7);
    clEnqueueSVMMap_ExpectAndReturn(commandQueuePool[0].get(), CL_TRUE, CL_MAP_READ|CL_MAP_WRITE, static_cast<void*>(vec.data()), vec.size()*sizeof(int), 0, NULL, NULL, CL_SUCCESS);
    TEST_ASSERT_EQUAL(commandQueuePool[0].enqueueMapSVM(vec, CL_TRUE, CL_MAP_READ|CL_MAP_WRITE), CL_SUCCESS);
#endif
}

// Run after other tests to clear the default state in the header
// using special unit test bypasses.
// We cannot remove the once_flag, so this is a hard fix
// but it means we won't hit cmock release callbacks at the end.
// This is a lot like tearDown but for the header default
// so we do not want to run it for every test.
// The alternative would be to manually modify the test runner
// but we avoid that for now.
void testCleanupHeaderState()
{
    clReleaseCommandQueue_ExpectAndReturn(make_command_queue(1), CL_SUCCESS);
    clReleaseContext_ExpectAndReturn(make_context(1), CL_SUCCESS);
    clReleaseDevice_ExpectAndReturn(make_device_id(1), CL_SUCCESS);

    cl::CommandQueue::unitTestClearDefault();
    cl::Context::unitTestClearDefault();
    cl::Device::unitTestClearDefault();
    cl::Platform::unitTestClearDefault();
}

// OpenCL 2.2 APIs:

static void CL_CALLBACK test_program_release_callback(
    cl_program,
    void*)
{
}

static cl_int clSetProgramReleaseCallback_set(
    cl_program program,
    void (CL_CALLBACK * pfn_notify)(cl_program program, void * user_data),
    void *user_data,
    int num_calls)
{
    (void) num_calls;

    TEST_ASSERT_EQUAL_PTR(make_program(0), program);
    TEST_ASSERT_EQUAL_PTR(pfn_notify, test_program_release_callback);

    return CL_SUCCESS;
}

void testSetProgramReleaseCallback()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 220
    cl_program program = make_program(0);
    int user_data = 0;

    clSetProgramReleaseCallback_StubWithCallback(clSetProgramReleaseCallback_set);
    clReleaseProgram_ExpectAndReturn(program, CL_SUCCESS);

    cl::Program prog(program);

    prog.setReleaseCallback(test_program_release_callback, &user_data);
#endif
}

void testSetProgramSpecializationConstantScalar()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 220
    cl_program program = make_program(0);
    int sc = 0;

    clSetProgramSpecializationConstant_ExpectAndReturn(program, 0, sizeof(sc), &sc, CL_SUCCESS);
    clReleaseProgram_ExpectAndReturn(program, CL_SUCCESS);

    cl::Program prog(program);

    prog.setSpecializationConstant(0, sc);
#endif
}

/// Stub for testing boolean specialization constants
static cl_int clSetProgramSpecializationConstant_testBool(
    cl_program program,
    cl_uint spec_id,
    size_t spec_size,
    const void* spec_value,
    int num_calls)
{
    (void) num_calls;

    TEST_ASSERT_EQUAL_PTR(make_program(0), program);
    TEST_ASSERT(spec_id == 0 || spec_id == 1);
    TEST_ASSERT_EQUAL(spec_size, 1);
    if (spec_id == 0)
    {
        const cl_uchar *uc_value = (const cl_uchar*)spec_value;
        TEST_ASSERT_EQUAL_HEX(uc_value[0], 0);
    }
    if (spec_id == 1)
    {
        const cl_uchar *uc_value = (const cl_uchar*)spec_value;
        TEST_ASSERT_EQUAL_HEX(uc_value[0], CL_UCHAR_MAX);
    }
    return CL_SUCCESS;
}

void testSetProgramSpecializationConstantBool()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 220
    // Spec constant "false" should turn into a call with size one and no bits set.
    // Spec constant "true" should turn into a call with size one and all bits set.
    cl_program program = make_program(0);
    bool scFalse = false;
    bool scTrue = true;

    clSetProgramSpecializationConstant_StubWithCallback(clSetProgramSpecializationConstant_testBool);

    clReleaseProgram_ExpectAndReturn(program, CL_SUCCESS);

    cl::Program prog(program);

    prog.setSpecializationConstant(0, scFalse);
    prog.setSpecializationConstant(1, scTrue);
#endif
}

void testSetProgramSpecializationConstantPointer()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 220
    cl_program program = make_program(0);
    int scArray[5];

    clSetProgramSpecializationConstant_ExpectAndReturn(program, 0, sizeof(scArray), &scArray, CL_SUCCESS);
    clReleaseProgram_ExpectAndReturn(program, CL_SUCCESS);

    cl::Program prog(program);

    prog.setSpecializationConstant(0, sizeof(scArray), scArray);
#endif
}

// OpenCL 3.0 and cl_khr_extended_versioning Queries

// Assumes the core enums, structures, and macros exactly match
// the extension enums, structures, and macros:

static_assert(CL_PLATFORM_NUMERIC_VERSION == CL_PLATFORM_NUMERIC_VERSION_KHR,
    "CL_PLATFORM_NUMERIC_VERSION mismatch");
static_assert(CL_PLATFORM_EXTENSIONS_WITH_VERSION == CL_PLATFORM_EXTENSIONS_WITH_VERSION_KHR,
    "CL_PLATFORM_EXTENSIONS_WITH_VERSION mismatch");

static_assert(CL_DEVICE_NUMERIC_VERSION == CL_DEVICE_NUMERIC_VERSION_KHR,
    "CL_DEVICE_NUMERIC_VERSION mismatch");
static_assert(CL_DEVICE_EXTENSIONS_WITH_VERSION == CL_DEVICE_EXTENSIONS_WITH_VERSION_KHR,
    "CL_DEVICE_EXTENSIONS_WITH_VERSION mismatch");
static_assert(CL_DEVICE_ILS_WITH_VERSION == CL_DEVICE_ILS_WITH_VERSION_KHR,
    "CL_DEVICE_ILS_WITH_VERSION mismatch");
static_assert(CL_DEVICE_BUILT_IN_KERNELS_WITH_VERSION == CL_DEVICE_BUILT_IN_KERNELS_WITH_VERSION_KHR,
    "CL_DEVICE_BUILT_IN_KERNELS_WITH_VERSION mismatch");

static_assert(sizeof(cl_name_version) == sizeof(cl_name_version_khr),
    "cl_name_version mismatch");

static_assert(CL_MAKE_VERSION(1, 2, 3) == CL_MAKE_VERSION_KHR(1, 2, 3),
    "CL_MAKE_VERSION mismatch");

static cl_int clGetPlatformInfo_extended_versioning(
    cl_platform_id id,
    cl_platform_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void)num_calls;
    switch (param_name) {
    case CL_PLATFORM_NUMERIC_VERSION:
    {
        if (param_value_size == sizeof(cl_version) && param_value) {
            *static_cast<cl_version*>(param_value) = CL_MAKE_VERSION(1, 2, 3);
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(cl_version);
        }
        return CL_SUCCESS;
    }
    case CL_PLATFORM_EXTENSIONS_WITH_VERSION:
    {
        static cl_name_version extension = {
            CL_MAKE_VERSION(10, 11, 12),
            "cl_dummy_extension",
        };
        if (param_value_size == sizeof(cl_name_version) && param_value) {
            *static_cast<cl_name_version*>(param_value) = extension;
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(extension);
        }
        return CL_SUCCESS;
    }
    default: break;
    }
    TEST_FAIL();
    return CL_INVALID_OPERATION;
}

void testPlatformExtendedVersioning_3_0()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 300
    cl::Platform p(make_platform_id(1));

    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_extended_versioning);

    cl_version platformVersion = p.getInfo<CL_PLATFORM_NUMERIC_VERSION>();
    TEST_ASSERT_EQUAL_HEX(platformVersion, CL_MAKE_VERSION(1, 2, 3));

    std::vector<cl_name_version> extensions = p.getInfo<CL_PLATFORM_EXTENSIONS_WITH_VERSION>();
    TEST_ASSERT_EQUAL(extensions.size(), 1);
    TEST_ASSERT_EQUAL_HEX(extensions[0].version, CL_MAKE_VERSION(10, 11, 12));
    TEST_ASSERT_EQUAL_STRING(extensions[0].name, "cl_dummy_extension");
#endif // CL_HPP_TARGET_OPENCL_VERSION >= 300
}

void testPlatformExtendedVersioning_KHR()
{
#if CL_HPP_TARGET_OPENCL_VERSION < 300
    cl::Platform p(make_platform_id(1));

    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_extended_versioning);

    cl_version_khr platformVersion = p.getInfo<CL_PLATFORM_NUMERIC_VERSION_KHR>();
    TEST_ASSERT_EQUAL_HEX(platformVersion, CL_MAKE_VERSION_KHR(1, 2, 3));

    std::vector<cl_name_version_khr> extensions = p.getInfo<CL_PLATFORM_EXTENSIONS_WITH_VERSION_KHR>();
    TEST_ASSERT_EQUAL(extensions.size(), 1);
    TEST_ASSERT_EQUAL_HEX(extensions[0].version, CL_MAKE_VERSION_KHR(10, 11, 12));
    TEST_ASSERT_EQUAL_STRING(extensions[0].name, "cl_dummy_extension");
#endif // CL_HPP_TARGET_OPENCL_VERSION < 300
}


// Note: This assumes the core enums, structures, and macros exactly match
// the extension enums, structures, and macros.

static cl_int clGetDeviceInfo_extended_versioning(
    cl_device_id id,
    cl_device_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void)num_calls;
    switch (param_name) {
    case CL_DEVICE_NUMERIC_VERSION:
    {
        if (param_value_size == sizeof(cl_version) && param_value) {
            *static_cast<cl_version*>(param_value) = CL_MAKE_VERSION(1, 2, 3);
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(cl_version);
        }
        return CL_SUCCESS;
    }
    case CL_DEVICE_OPENCL_C_NUMERIC_VERSION_KHR:
    {
        if (param_value_size == sizeof(cl_version_khr) && param_value) {
            *static_cast<cl_version_khr*>(param_value) = CL_MAKE_VERSION_KHR(4, 5, 6);
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(cl_version_khr);
        }
        return CL_SUCCESS;
    }
    case CL_DEVICE_EXTENSIONS_WITH_VERSION:
    {
        static cl_name_version extension = {
            CL_MAKE_VERSION(10, 11, 12),
            "cl_dummy_extension",
        };
        if (param_value_size == sizeof(cl_name_version) && param_value) {
            *static_cast<cl_name_version*>(param_value) = extension;
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(extension);
        }
        return CL_SUCCESS;
    }
    case CL_DEVICE_ILS_WITH_VERSION:
    {
        static cl_name_version il = {
            CL_MAKE_VERSION(20, 21, 22),
            "DUMMY_IR",
        };
        if (param_value_size == sizeof(cl_name_version) && param_value) {
            *static_cast<cl_name_version*>(param_value) = il;
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(il);
        }
        return CL_SUCCESS;
    }
    case CL_DEVICE_BUILT_IN_KERNELS_WITH_VERSION:
    {
        // Test no built-in kernels:
        if (param_value_size_ret) {
            *param_value_size_ret = 0;
        }
        return CL_SUCCESS;
    }
    case CL_DEVICE_OPENCL_C_ALL_VERSIONS:
    {
        static cl_name_version opencl_c = {
            CL_MAKE_VERSION(30, 31, 32),
            "OpenCL C",
        };
        if (param_value_size == sizeof(cl_name_version) && param_value) {
            *static_cast<cl_name_version*>(param_value) = opencl_c;
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(opencl_c);
        }
        return CL_SUCCESS;
    }
    case CL_DEVICE_OPENCL_C_FEATURES:
    {
        static cl_name_version opencl_c_features[] = {
            {
                CL_MAKE_VERSION(40, 41, 42),
                "__opencl_c_feature",
            },
            {
                CL_MAKE_VERSION(40, 43, 44),
                "__opencl_c_fancy_feature",
            },
        };
        if (param_value_size == sizeof(opencl_c_features) && param_value) {
            cl_name_version* feature = static_cast<cl_name_version*>(param_value);
            const int numFeatures = ARRAY_SIZE(opencl_c_features);
            for (int i = 0; i < numFeatures; i++) {
                feature[i] = opencl_c_features[i];
            }
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(opencl_c_features);
        }
        return CL_SUCCESS;
    }
    default: break;
    }
    TEST_FAIL();
    return CL_INVALID_OPERATION;
}

void testDeviceExtendedVersioning_3_0()
{
#if CL_HPP_TARGET_OPENCL_VERSION >= 300
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_3_0);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    cl::Device d0(make_device_id(0));

    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_extended_versioning);

    cl_version deviceVersion = d0.getInfo<CL_DEVICE_NUMERIC_VERSION>();
    TEST_ASSERT_EQUAL_HEX(deviceVersion, CL_MAKE_VERSION(1, 2, 3));

    std::vector<cl_name_version> extensions = d0.getInfo<CL_DEVICE_EXTENSIONS_WITH_VERSION>();
    TEST_ASSERT_EQUAL(extensions.size(), 1);
    TEST_ASSERT_EQUAL_HEX(extensions[0].version, CL_MAKE_VERSION(10, 11, 12));
    TEST_ASSERT_EQUAL_STRING(extensions[0].name, "cl_dummy_extension");

    std::vector<cl_name_version> ils = d0.getInfo<CL_DEVICE_ILS_WITH_VERSION>();
    TEST_ASSERT_EQUAL(ils.size(), 1);
    TEST_ASSERT_EQUAL_HEX(ils[0].version, CL_MAKE_VERSION(20, 21, 22));
    TEST_ASSERT_EQUAL_STRING(ils[0].name, "DUMMY_IR");

    std::vector<cl_name_version> opencl_c = d0.getInfo<CL_DEVICE_OPENCL_C_ALL_VERSIONS>();
    TEST_ASSERT_EQUAL(opencl_c.size(), 1);
    TEST_ASSERT_EQUAL_HEX(opencl_c[0].version, CL_MAKE_VERSION(30, 31, 32));
    TEST_ASSERT_EQUAL_STRING(opencl_c[0].name, "OpenCL C");

    std::vector<cl_name_version> opencl_c_features = d0.getInfo<CL_DEVICE_OPENCL_C_FEATURES>();
    TEST_ASSERT_EQUAL(opencl_c_features.size(), 2);
    TEST_ASSERT_EQUAL_HEX(opencl_c_features[0].version, CL_MAKE_VERSION(40, 41, 42));
    TEST_ASSERT_EQUAL_STRING(opencl_c_features[0].name, "__opencl_c_feature");
    TEST_ASSERT_EQUAL_HEX(opencl_c_features[1].version, CL_MAKE_VERSION(40, 43, 44));
    TEST_ASSERT_EQUAL_STRING(opencl_c_features[1].name, "__opencl_c_fancy_feature");

    std::vector<cl_name_version> builtInKernels = d0.getInfo<CL_DEVICE_BUILT_IN_KERNELS_WITH_VERSION>();
    TEST_ASSERT_EQUAL(builtInKernels.size(), 0);
#endif // CL_HPP_TARGET_OPENCL_VERSION >= 300
}

void testDeviceExtendedVersioning_KHR()
{
#if CL_HPP_TARGET_OPENCL_VERSION < 300
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    cl::Device d0(make_device_id(0));

    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_extended_versioning);

    cl_version_khr deviceVersion = d0.getInfo<CL_DEVICE_NUMERIC_VERSION_KHR>();
    TEST_ASSERT_EQUAL_HEX(deviceVersion, CL_MAKE_VERSION_KHR(1, 2, 3));

    cl_version_khr cVersion = d0.getInfo<CL_DEVICE_OPENCL_C_NUMERIC_VERSION_KHR>();
    TEST_ASSERT_EQUAL_HEX(cVersion, CL_MAKE_VERSION_KHR(4, 5, 6));

    std::vector<cl_name_version_khr> extensions = d0.getInfo<CL_DEVICE_EXTENSIONS_WITH_VERSION_KHR>();
    TEST_ASSERT_EQUAL(extensions.size(), 1);
    TEST_ASSERT_EQUAL_HEX(extensions[0].version, CL_MAKE_VERSION_KHR(10, 11, 12));
    TEST_ASSERT_EQUAL_STRING(extensions[0].name, "cl_dummy_extension");

    std::vector<cl_name_version_khr> ils = d0.getInfo<CL_DEVICE_ILS_WITH_VERSION_KHR>();
    TEST_ASSERT_EQUAL(ils.size(), 1);
    TEST_ASSERT_EQUAL_HEX(ils[0].version, CL_MAKE_VERSION_KHR(20, 21, 22));
    TEST_ASSERT_EQUAL_STRING(ils[0].name, "DUMMY_IR");

    std::vector<cl_name_version_khr> builtInKernels = d0.getInfo<CL_DEVICE_BUILT_IN_KERNELS_WITH_VERSION_KHR>();
    TEST_ASSERT_EQUAL(builtInKernels.size(), 0);
#endif // CL_HPP_TARGET_OPENCL_VERSION < 300
}

static cl_int clGetDeviceInfo_uuid_pci_bus_info(
    cl_device_id id,
    cl_device_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret,
    int num_calls)
{
    (void)num_calls;
    switch (param_name) {
#if defined(cl_khr_device_uuid)
    case CL_DEVICE_UUID_KHR:
    case CL_DRIVER_UUID_KHR:
    {
        if (param_value_size == CL_UUID_SIZE_KHR && param_value) {
            cl_uchar* pUUID = static_cast<cl_uchar*>(param_value);
            cl_uchar start =
                (param_name == CL_DEVICE_UUID_KHR) ? 1 :
                (param_name == CL_DRIVER_UUID_KHR) ? 2 :
                0;
            for (int i = 0; i < CL_UUID_SIZE_KHR; i++) {
                pUUID[i] = i + start;
            }
        }
        if (param_value_size_ret) {
            *param_value_size_ret = CL_UUID_SIZE_KHR;
        }
        return CL_SUCCESS;
    }
    case CL_DEVICE_LUID_VALID_KHR:
    {
        if (param_value_size == sizeof(cl_bool) && param_value) {
            *static_cast<cl_bool*>(param_value) = CL_TRUE;
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(cl_bool);
        }
        return CL_SUCCESS;
    }
    case CL_DEVICE_LUID_KHR:
    {
        if (param_value_size == CL_LUID_SIZE_KHR && param_value) {
            cl_uchar* pLUID = static_cast<cl_uchar*>(param_value);
            cl_uchar start = 3;
            for (int i = 0; i < CL_LUID_SIZE_KHR; i++) {
                pLUID[i] = i + start;
            }
        }
        if (param_value_size_ret) {
            *param_value_size_ret = CL_LUID_SIZE_KHR;
        }
        return CL_SUCCESS;
    }
    case CL_DEVICE_NODE_MASK_KHR:
    {
        if (param_value_size == sizeof(cl_uint) && param_value) {
            *static_cast<cl_uint*>(param_value) = 0xA5A5;
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(cl_uint);
        }
        return CL_SUCCESS;
    }
#endif
#if defined(cl_khr_pci_bus_info)
    case CL_DEVICE_PCI_BUS_INFO_KHR:
    {
        if (param_value_size == sizeof(cl_device_pci_bus_info_khr) && param_value) {
            cl_device_pci_bus_info_khr* pInfo = static_cast<cl_device_pci_bus_info_khr*>(param_value);
            pInfo->pci_domain = 0x11;
            pInfo->pci_bus = 0x22;
            pInfo->pci_device = 0x33;
            pInfo->pci_function = 0x44;
        }
        if (param_value_size_ret) {
            *param_value_size_ret = sizeof(cl_device_pci_bus_info_khr);
        }
        return CL_SUCCESS;
    }
#endif
    default: break;
    }
    TEST_FAIL();
    return CL_INVALID_OPERATION;
}

void testDeviceUUID_KHR()
{
#if defined(cl_khr_device_uuid)
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    cl::Device d0(make_device_id(0));

    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_uuid_pci_bus_info);

    std::array<cl_uchar, CL_UUID_SIZE_KHR> dev_uuid = d0.getInfo<CL_DEVICE_UUID_KHR>();
    for (int i = 0; i < CL_UUID_SIZE_KHR; i++) {
        TEST_ASSERT_EQUAL_UINT8(i + 1, dev_uuid[i]);
    }
    std::array<cl_uchar, CL_UUID_SIZE_KHR> drv_uuid = d0.getInfo<CL_DRIVER_UUID_KHR>();
    for (int i = 0; i < CL_UUID_SIZE_KHR; i++) {
        TEST_ASSERT_EQUAL_UINT8(i + 2, drv_uuid[i]);
    }

    cl_bool valid = d0.getInfo<CL_DEVICE_LUID_VALID_KHR>();
    TEST_ASSERT_EQUAL(CL_TRUE, valid);
    std::array<cl_uchar, CL_LUID_SIZE_KHR> luid = d0.getInfo<CL_DEVICE_LUID_KHR>();
    for (int i = 0; i < CL_LUID_SIZE_KHR; i++) {
        TEST_ASSERT_EQUAL_UINT8(i + 3, luid[i]);
    }

    cl_uint nodeMask = d0.getInfo<CL_DEVICE_NODE_MASK_KHR>();
    TEST_ASSERT_EQUAL(0xA5A5, nodeMask);
#endif
}

void testDevicePCIBusInfo_KHR()
{
#if defined(cl_khr_pci_bus_info)
    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_platform);
    clGetPlatformInfo_StubWithCallback(clGetPlatformInfo_version_2_0);
    clReleaseDevice_ExpectAndReturn(make_device_id(0), CL_SUCCESS);

    cl::Device d0(make_device_id(0));

    clGetDeviceInfo_StubWithCallback(clGetDeviceInfo_uuid_pci_bus_info);

    cl_device_pci_bus_info_khr info = d0.getInfo<CL_DEVICE_PCI_BUS_INFO_KHR>();
    TEST_ASSERT_EQUAL_HEX(0x11, info.pci_domain);
    TEST_ASSERT_EQUAL_HEX(0x22, info.pci_bus);
    TEST_ASSERT_EQUAL_HEX(0x33, info.pci_device);
    TEST_ASSERT_EQUAL_HEX(0x44, info.pci_function);
#endif
}


} // extern "C"
