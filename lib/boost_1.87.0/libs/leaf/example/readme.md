# Example Programs Using LEAF to Handle Errors

* [print_file](./print_file): This directory contains several versions of a trivial program which takes a file name on the command line and prints it. Each version uses a different error handling implementaiton.

* [try_capture_all_result.cpp](https://github.com/boostorg/leaf/blob/master/example/try_capture_all_result.cpp?ts=4): Shows how to transport error objects between threads in a `leaf::result<T>` object without using exception handling.
* [try_capture_all_exceptions.cpp](https://github.com/boostorg/leaf/blob/master/example/try_capture_all_exceptions.cpp?ts=4): Shows how to transport error objects between threads in a `leaf::result<T>` object using exception handling.
* [lua_callback_result.cpp](https://github.com/boostorg/leaf/blob/master/example/lua_callback_result.cpp?ts=4): Transporting arbitrary error objects through an uncooperative C API.
* [lua_callback_exceptions.cpp](https://github.com/boostorg/leaf/blob/master/example/lua_callback_exceptions.cpp?ts=4): Transporting arbitrary error objects through an uncooperative API using exceptions.
* [exception_to_result.cpp](https://github.com/boostorg/leaf/blob/master/example/exception_to_result.cpp?ts=4): Demonstrates how to transport exceptions through a `noexcept` layer in the program.
* [exception_error_log.cpp](https://github.com/boostorg/leaf/blob/master/example/error_log.cpp?ts=4): Using `accumulate` to produce an error log.
* [exception_error_trace.cpp](https://github.com/boostorg/leaf/blob/master/example/error_trace.cpp?ts=4): Same as above, but the log is recorded in a `std::deque` rather than just printed.
* [print_half.cpp](https://github.com/boostorg/leaf/blob/master/example/print_half.cpp?ts=4): This is a Boost Outcome example adapted to LEAF, demonstrating the use of `try_handle_some` to handle some errors, forwarding any other error to the caller.
