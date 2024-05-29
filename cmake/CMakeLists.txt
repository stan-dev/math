set(BUILD_TESTING OFF)
set(BUILD_TESTING_STATIC OFF)
set(BUILD_TESTING_SHARED OFF)

# Externally provided libraries
FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG main)
FetchContent_Declare(googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG main) # need master for benchmark::benchmark

FetchContent_MakeAvailable(googletest googlebenchmark)


if(STAN_USE_SYSTEM_TBB)
  find_package(TBB REQUIRED)
else()
  FetchContent_Declare(
    tbb
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
    GIT_TAG        v2021.7.0  # adjust this to the version you need
  )
  FetchContent_MakeAvailable(tbb)
endif()

if(STAN_USE_SYSTEM_EIGEN)
  find_package(Eigen3 REQUIRED)
else()
  set(EIGEN_BUILD_DOC OFF)
  # note: To disable eigen tests,
  # you should put this code in a add_subdirectory to avoid to change
  # BUILD_TESTING for your own project too since variables are directory
  # scoped
  set(EIGEN_BUILD_BTL OFF)
  set(EIGEN_BUILD_TESTING OFF)
  set(EIGEN_BUILD_PKGCONFIG OFF)
  set(EIGEN_LEAVE_TEST_IN_ALL_TARGET OFF)
  FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0)
  FetchContent_MakeAvailable(Eigen)
  #TODO: Need to just put dependency file in subdirectory
  # since macros are directory scoped
  set(BUILD_TESTING ON)
  set(EIGEN_BUILD_TESTING ON)
endif()

if (STAN_USE_SYSTEM_SUNDIALS)
  find_package(SUNDIALS REQUIRED)
else()
  FetchContent_Declare(
    sundials
    GIT_REPOSITORY https://github.com/LLNL/sundials
    GIT_TAG        v6.1.1
      # adjust this to the version you need
  )
  FetchContent_GetProperties(sundials)
  if(NOT sundials_POPULATED)
    FetchContent_Populate(sundials)
    add_subdirectory(${sundials_SOURCE_DIR} ${sundials_BINARY_DIR})
  endif()
endif()

if (STAN_USE_SYSTEM_BOOST)
  find_package(Boost REQUIRED)
else()
  set(Boost_USE_STATIC_LIBS ON)
  set(Boost_USE_STATIC_RUNTIME ON)
  set(BOOST_ENABLE_CMAKE ON)
  FetchContent_Declare(
    boost
    URL https://boostorg.jfrog.io/artifactory/main/release/1.85.0/source/boost_1_85_0.tar.gz
  )
  FetchContent_GetProperties(boost)
  if(NOT boost_POPULATED)
    FetchContent_Populate(boost)
    set(boost_INCLUDE_DIR ${boost_SOURCE_DIR})
  endif()
endif()
