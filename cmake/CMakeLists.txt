# Share dependency versions with Make; cache paths may be overridden.
file(STRINGS "${PROJECT_SOURCE_DIR}/make/libraries" library_locations
  REGEX "^(BOOST|EIGEN|OPENCL|TBB|SUNDIALS|BENCHMARK) ")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/make/libraries")
foreach(line IN LISTS library_locations)
  if(line MATCHES "^([A-Z]+) .*lib/([^ ]+)")
    set(dependency "${CMAKE_MATCH_1}")
    set(vendored_root "${PROJECT_SOURCE_DIR}/lib/${CMAKE_MATCH_2}")
    set(STAN_${dependency}_ROOT "" CACHE PATH "Override ${dependency} source directory (empty: vendored)")
    if(NOT STAN_${dependency}_ROOT)
      set(STAN_${dependency}_ROOT "${vendored_root}")
    endif()
    set(STAN_${dependency}_ROOT "${STAN_${dependency}_ROOT}" PARENT_SCOPE)
  endif()
endforeach()
option(STAN_USE_SYSTEM_TBB "Use an installed TBB library" OFF)
option(STAN_USE_SYSTEM_EIGEN "Use an installed Eigen library" OFF)
option(STAN_USE_SYSTEM_BOOST "Use installed Boost headers" OFF)
option(STAN_USE_SYSTEM_SUNDIALS "Use installed SUNDIALS libraries" OFF)
find_package(Threads REQUIRED)
target_link_libraries(stanmath INTERFACE Threads::Threads)

if(STAN_USE_SYSTEM_EIGEN)
  find_package(Eigen3 REQUIRED)
  target_link_libraries(stanmath INTERFACE Eigen3::Eigen)
else()
  target_include_directories(stanmath SYSTEM INTERFACE "${STAN_EIGEN_ROOT}")
endif()
if(STAN_USE_SYSTEM_BOOST)
  find_package(Boost REQUIRED)
  target_link_libraries(stanmath INTERFACE Boost::boost)
else()
  target_include_directories(stanmath SYSTEM INTERFACE "${STAN_BOOST_ROOT}")
endif()

if(STAN_USE_SYSTEM_TBB)
  find_package(TBB REQUIRED COMPONENTS tbb)
else()
  # Vendored TBB predates oneTBB's CMake build. One nested Make job keeps
  # dependency compilation within the outer build's parallelism limit.
  include(ExternalProject)
  find_program(STAN_TBB_MAKE NAMES gmake make)
  if(NOT STAN_TBB_MAKE)
    message(FATAL_ERROR "Building vendored TBB requires GNU Make")
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(tbb_compiler clang)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(tbb_compiler gcc)
  else()
    message(FATAL_ERROR "Use STAN_USE_SYSTEM_TBB with this compiler")
  endif()
  if(APPLE)
    set(tbb_library "${CMAKE_CURRENT_BINARY_DIR}/tbb/libtbb.dylib")
  elseif(WIN32)
    message(FATAL_ERROR "Use STAN_USE_SYSTEM_TBB on Windows")
  else()
    set(tbb_library "${CMAKE_CURRENT_BINARY_DIR}/tbb/libtbb.so.2")
  endif()
  ExternalProject_Add(stan_tbb_build
    SOURCE_DIR "${STAN_TBB_ROOT}"
    BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/tbb"
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${STAN_TBB_MAKE} -j1 -r -f "${STAN_TBB_ROOT}/build/Makefile.tbb"
      "tbb_root=${STAN_TBB_ROOT}" compiler=${tbb_compiler} cfg=release stdver=c++17
      "CXX=${CMAKE_CXX_COMPILER}" "CC=${CMAKE_C_COMPILER}"
      "CXXFLAGS=-Wno-deprecated-declarations -Wno-deprecated-copy"
    BUILD_BYPRODUCTS "${tbb_library}"
    BUILD_ALWAYS TRUE
    INSTALL_COMMAND "")
  add_library(TBB::tbb SHARED IMPORTED GLOBAL)
  set_target_properties(TBB::tbb PROPERTIES
    IMPORTED_LOCATION "${tbb_library}"
    INTERFACE_INCLUDE_DIRECTORIES "${STAN_TBB_ROOT}/include")
  add_dependencies(TBB::tbb stan_tbb_build)
endif()
target_link_libraries(stanmath INTERFACE TBB::tbb)

if(STAN_USE_SYSTEM_SUNDIALS)
  find_package(SUNDIALS REQUIRED COMPONENTS cvodes idas kinsol nvecserial)
  target_link_libraries(stanmath INTERFACE
    SUNDIALS::cvodes SUNDIALS::idas SUNDIALS::kinsol SUNDIALS::nvecserial)
else()
  # The pruned vendored distribution has no upstream CMakeLists.txt.
  # Match make/libraries, including Stan's printf override and common modules.
  file(GLOB sundials_sources CONFIGURE_DEPENDS
    "${STAN_SUNDIALS_ROOT}/src/cvodes/*.c"
    "${STAN_SUNDIALS_ROOT}/src/idas/*.c"
    "${STAN_SUNDIALS_ROOT}/src/kinsol/*.c"
    "${STAN_SUNDIALS_ROOT}/src/sundials/*.c"
    "${STAN_SUNDIALS_ROOT}/src/sunmatrix/band/*.c"
    "${STAN_SUNDIALS_ROOT}/src/sunmatrix/dense/*.c"
    "${STAN_SUNDIALS_ROOT}/src/sunlinsol/band/*.c"
    "${STAN_SUNDIALS_ROOT}/src/sunlinsol/dense/*.c"
    "${STAN_SUNDIALS_ROOT}/src/sunnonlinsol/newton/*.c"
    "${STAN_SUNDIALS_ROOT}/src/sunnonlinsol/fixedpoint/*.c"
    "${STAN_SUNDIALS_ROOT}/src/nvector/serial/nvector_serial.c")
  list(FILTER sundials_sources EXCLUDE REGEX "/f[^/]*\\.c$|/sundials_profiler\\.c$")
  add_library(stan_sundials STATIC ${sundials_sources})
  if(BUILD_TESTING)
    target_compile_options(stan_sundials PRIVATE ${STAN_TEST_OPTIMIZATION_FLAGS})
  endif()
  target_include_directories(stan_sundials SYSTEM PUBLIC "${STAN_SUNDIALS_ROOT}/include")
  target_include_directories(stan_sundials PRIVATE "${STAN_SUNDIALS_ROOT}/src/sundials")
  target_compile_definitions(stan_sundials PRIVATE NO_FPRINTF_OUTPUT)
  if(MSVC)
    target_compile_options(stan_sundials PRIVATE
      "/FI${STAN_SUNDIALS_ROOT}/include/stan_sundials_printf_override.hpp")
  else()
    target_compile_options(stan_sundials PRIVATE
      -include "${STAN_SUNDIALS_ROOT}/include/stan_sundials_printf_override.hpp")
    target_link_libraries(stan_sundials PUBLIC m)
  endif()
  target_link_libraries(stanmath INTERFACE stan_sundials)
endif()
if(BUILD_TESTING)
  set(BUILD_GMOCK OFF CACHE BOOL "Build Google Mock")
  set(INSTALL_GTEST OFF CACHE BOOL "Install Google Test")
  add_subdirectory("${STAN_BENCHMARK_ROOT}/googletest"
    "${CMAKE_CURRENT_BINARY_DIR}/googletest" EXCLUDE_FROM_ALL)
endif()
