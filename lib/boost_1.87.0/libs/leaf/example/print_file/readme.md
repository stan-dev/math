# Print File Example

This directory contains several versions of a trivial program which takes a file name on the command line and prints it. Each version uses a different error handling implementaiton.

* [print_file_leaf_result.cpp](./print_file_leaf_result.cpp) reports errors with
  `leaf::result<T>`, using an error code `enum` for classification of failures.

* [print_file_system_result.cpp](./print_file_system_result.cpp) is the same as
  above, but using `boost::system::result<T>` instead of `leaf::result<T>`.
  This demonstrates the ability of LEAF to transport arbitrary error objects using an
  external result type, rather than `boost::leaf::result<T>`.

* [print_file_exceptions.cpp](./print_file_exceptions.cpp) throws on error, using an error code
  `enum` for classification of failures.
