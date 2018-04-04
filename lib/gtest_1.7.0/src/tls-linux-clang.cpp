// needed on linux when using clang for thread_local storage to work
// see: https://stackoverflow.com/questions/29322666/undefined-reference-to-cxa-thread-atexitcxxabi-when-compiling-with-libc#30437761
extern "C" int __cxa_thread_atexit(void (*func)(), void *obj,
                                   void *dso_symbol) {
  int __cxa_thread_atexit_impl(void (*)(), void *, void *);
  return __cxa_thread_atexit_impl(func, obj, dso_symbol);
}
