# Boost.MySQL

Branch | Windows/Linux Build | OSX build | Coverage | Documentation
-------|---------------------|-----------|--------- | -------------
[`master`](https://github.com/boostorg/mysql/tree/master)   | [![Build Status](https://drone.cpp.al/api/badges/boostorg/mysql/status.svg)](https://drone.cpp.al/boostorg/mysql)                        | [![Build Status](https://github.com/boostorg/mysql/actions/workflows/build-code.yml/badge.svg)](https://github.com/boostorg/mysql)                | [![codecov](https://codecov.io/gh/boostorg/mysql/branch/master/graph/badge.svg)](https://codecov.io/gh/boostorg/mysql/branch/master)   | [Docs for master](https://www.boost.org/doc/libs/master/libs/mysql/doc/html/index.html)
[`develop`](https://github.com/boostorg/mysql/tree/develop) | [![Build Status](https://drone.cpp.al/api/badges/boostorg/mysql/status.svg?ref=refs/heads/develop)](https://drone.cpp.al/boostorg/mysql) | [![Build Status](https://github.com/boostorg/mysql/actions/workflows/build-code.yml/badge.svg?branch=develop)](https://github.com/boostorg/mysql) | [![codecov](https://codecov.io/gh/boostorg/mysql/branch/develop/graph/badge.svg)](https://codecov.io/gh/boostorg/mysql/branch/develop) | [Docs for develop](https://www.boost.org/doc/libs/develop/libs/mysql/doc/html/index.html)

Boost.MySQL is a C++11 client for MySQL and MariaDB database servers, based on Boost.Asio.
Boost.MySQL is part of Boost.

Boost.MySQL gets its first stable release with Boost 1.82 on the 12th of April, 2023.

## Why another MySQL C++ client?

- It is fully compatible with Boost.Asio and integrates well with any other
  library in the Boost.Asio ecosystem (like Boost.Beast).
- It supports Boost.Asio's universal asynchronous model, which means you can
  go asyncrhonous using callbacks, futures or coroutines (including C++20 coroutines).
- It is written in C++11 and takes advantage of it.
- It is header only.

## Using the library

To use this library, you need:

- A C++11 capable compiler.
- Boost 1.82 or higher (Boost.MySQL doesn't work with standalone Asio).
- OpenSSL.

The library is header-only, but it depends on other Boost header-only libraries and on OpenSSL.
To use the library, install Boost the way you would normally do (e.g. via `b2 install`), and create
a `CMakeLists.txt` like this (replace `main` by your executable name and `main.cpp` by your list of source files):

```cmake
project(boost_mysql_example LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS headers)
find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE Boost::headers Threads::Threads OpenSSL::Crypto OpenSSL::SSL)
```

## Tested with

Boost.MySQL has been tested with the following compilers:
- gcc 5 to 11.
- clang 3.6 to 14.
- msvc 14.1, 14.2 and 14.3.

And with the following databases:
- MySQL v5.7.41.
- MySQL v8.0.31.
- MariaDB v10.11.2.

## Upgrading from 0.2.x

If you were using 0.2.x, you can find [upgrade instructions here](doc/upgrade_1_82.md). 

## Features

Currently implemented:
- Text queries (execution of text SQL queries and data retrieval).
  MySQL refers to this as the "text protocol", as all information is passed using text
  (as opposed to prepared statements, see below).
- Prepared statements. MySQL refers to this as the "binary protocol", as the result
  of executing a prepared statement is sent in binary format rather than in text.
- Authentication methods (authentication plugins): mysql_native_password and
  caching_sha2_password. These are the default methods in MySQL 5 and MySQL 8,
  respectively.
- Encrypted connections (TLS).
- TCP and UNIX socket connections. The implementation is based on Boost.Asio
  SyncStream and AsyncStream concepts, so it is generic and can be used with
  any stream that fulfills these concept's requirements. There are user-friendly
  typedefs and regression tests for TCP and UNIX socket streams.

Yet to be done (but it is on our list - PRs welcome):
- Further authentication methods: sha256_password
- Multi-resultset: being able to specify several semicolon-separated queries.
