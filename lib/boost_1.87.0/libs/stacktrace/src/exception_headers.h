// Copyright Antony Polukhin, 2023-2024.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <stddef.h>

#if defined(__x86_64__) || defined(_M_X64) || defined(__MINGW32__)
#define BOOST_STACKTRACE_ALWAYS_STORE_IN_PADDING 1
#else
#define BOOST_STACKTRACE_ALWAYS_STORE_IN_PADDING 0
#endif


extern "C" {

// Developer note: helper to experiment with layouts of different
// exception headers https://godbolt.org/z/rrcdPbh1P

// https://github.com/llvm/llvm-project/blob/b3dd14ce07f2750ae1068fe62abbf2f3bd2cade8/libcxxabi/src/cxa_exception.h
struct cxa_exception_begin_llvm {
    const char* reserve;
    size_t referenceCount;
};

static cxa_exception_begin_llvm* exception_begin_llvm_ptr(void* ptr) {
  size_t kExceptionBeginOffset = (
    sizeof(void*) == 8 ? 128 : 80
  );
  return (cxa_exception_begin_llvm*)((char*)ptr - kExceptionBeginOffset);
}

// https://github.com/gcc-mirror/gcc/blob/5d2a360f0a541646abb11efdbabc33c6a04de7ee/libstdc%2B%2B-v3/libsupc%2B%2B/unwind-cxx.h#L100
struct cxa_exception_begin_gcc {
    size_t referenceCount;
    const char* reserve;
};

static cxa_exception_begin_gcc* exception_begin_gcc_ptr(void* ptr) {
  size_t kExceptionBeginOffset = (
    sizeof(void*) == 8 ? 128 : 96
  );
  return (cxa_exception_begin_gcc*)((char*)ptr - kExceptionBeginOffset);
}

static void* get_current_exception_raw_ptr(void* exc_ptr) {
  // https://github.com/gcc-mirror/gcc/blob/16e2427f50c208dfe07d07f18009969502c25dc8/libstdc%2B%2B-v3/libsupc%2B%2B/eh_ptr.cc#L147
  return *static_cast<void**>(exc_ptr);
}

}  // extern "C"


