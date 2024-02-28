# Copyright 2020, 2021 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.11...3.20)

project(boost_locale VERSION "${BOOST_SUPERPROJECT_VERSION}" LANGUAGES CXX)

# Using glob here is ok as it is only for headers
file(GLOB_RECURSE headers include/*.hpp)

add_library(boost_locale
  src/boost/locale/encoding/codepage.cpp
  src/boost/locale/encoding/iconv_converter.hpp
  src/boost/locale/encoding/uconv_converter.hpp
  src/boost/locale/encoding/wconv_converter.hpp
  src/boost/locale/shared/date_time.cpp
  src/boost/locale/shared/format.cpp
  src/boost/locale/shared/formatting.cpp
  src/boost/locale/shared/generator.cpp
  src/boost/locale/shared/iconv_codecvt.cpp
  src/boost/locale/shared/iconv_codecvt.hpp
  src/boost/locale/shared/ids.cpp
  src/boost/locale/shared/localization_backend.cpp
  src/boost/locale/shared/message.cpp
  src/boost/locale/shared/mo_lambda.cpp
  src/boost/locale/util/codecvt_converter.cpp
  src/boost/locale/util/default_locale.cpp
  src/boost/locale/util/encoding.cpp
  src/boost/locale/util/encoding.hpp
  src/boost/locale/util/foreach_char.hpp
  src/boost/locale/util/info.cpp
  src/boost/locale/util/locale_data.cpp
  src/boost/locale/util/make_std_unique.hpp
  src/boost/locale/util/numeric.hpp
  src/boost/locale/util/timezone.hpp
  ${headers}
)

add_library(Boost::locale ALIAS boost_locale)

target_include_directories(boost_locale PUBLIC include PRIVATE src)
target_compile_features(boost_locale PUBLIC cxx_std_11)

target_link_libraries(boost_locale
  PUBLIC
    Boost::assert
    Boost::config
    Boost::core
    Boost::iterator
    Boost::utility
  PRIVATE
    Boost::predef
    Boost::thread
)

target_compile_definitions(boost_locale
  PUBLIC
    BOOST_LOCALE_NO_LIB
  PRIVATE
    BOOST_LOCALE_SOURCE
    _CRT_SECURE_NO_WARNINGS
    _SCL_SECURE_NO_WARNINGS
)

if(BUILD_SHARED_LIBS)
  target_compile_definitions(boost_locale PUBLIC BOOST_LOCALE_DYN_LINK)
else()
  target_compile_definitions(boost_locale PUBLIC BOOST_LOCALE_STATIC_LINK)
endif()

# Build Options

find_package(ICU COMPONENTS data i18n uc QUIET)
option(BOOST_LOCALE_ENABLE_ICU "Boost.Locale: enable ICU backend" ${ICU_FOUND})

# Workaround for an issue e.g. on FreeBSD where including ICU may change the used libiconv
if(UNIX AND NOT DEFINED Iconv_INCLUDE_DIR AND BOOST_LOCALE_ENABLE_ICU AND EXISTS "${ICU_INCLUDE_DIRS}/iconv.h")
  set(Iconv_INCLUDE_DIR "${ICU_INCLUDE_DIRS}" CACHE PATH "Path to iconv headers")
  find_package(Iconv QUIET)
  if(NOT Iconv_FOUND)
    # Unset and retry below
    unset(Iconv_INCLUDE_DIR CACHE)
    unset(Iconv_LIBRARY CACHE)
    unset(Iconv_IS_BUILT_IN CACHE)
    unset(Iconv_IS_BUILT_IN)
  endif()
endif()

find_package(Iconv QUIET)
option(BOOST_LOCALE_ENABLE_ICONV "Boost.Locale: enable iconv backend" ${Iconv_FOUND})

if(CMAKE_SYSTEM_NAME MATCHES "Linux|Darwin")
  set(_default_posix ON)
elseif(CMAKE_SYSTEM_NAME MATCHES "FreeBSD")
  include (CheckIncludeFileCXX)
  check_include_file_cxx ("xlocale.h" HAVE_XLOCALE_H)
  set(_default_posix ${HAVE_XLOCALE_H})
else()
  set(_default_posix OFF)
endif()

option(BOOST_LOCALE_ENABLE_POSIX "Boost.Locale: enable POSIX backend" ${_default_posix})
unset(_default_posix)

option(BOOST_LOCALE_ENABLE_STD "Boost.Locale: enable std::locale backend" ON)
option(BOOST_LOCALE_ENABLE_WINAPI "Boost.Locale: enable Windows API backend" ${WIN32})


if(BOOST_LOCALE_ENABLE_ICONV)

  find_package(Iconv REQUIRED)

  target_compile_definitions(boost_locale PRIVATE BOOST_LOCALE_WITH_ICONV=1)
  target_link_libraries(boost_locale PRIVATE Iconv::Iconv)

  target_sources(boost_locale PRIVATE
    src/boost/locale/util/iconv.hpp
  )
endif()

if(BOOST_LOCALE_ENABLE_ICU)

  find_package(ICU COMPONENTS data i18n uc REQUIRED)

  target_compile_definitions(boost_locale PRIVATE BOOST_LOCALE_WITH_ICU=1)
  target_link_libraries(boost_locale PRIVATE ICU::data ICU::i18n ICU::uc)

  target_sources(boost_locale PRIVATE
    src/boost/locale/icu/all_generator.hpp
    src/boost/locale/icu/boundary.cpp
    src/boost/locale/icu/cdata.hpp
    src/boost/locale/icu/codecvt.cpp
    src/boost/locale/icu/codecvt.hpp
    src/boost/locale/icu/collator.cpp
    src/boost/locale/icu/conversion.cpp
    src/boost/locale/icu/date_time.cpp
    src/boost/locale/icu/formatter.cpp
    src/boost/locale/icu/formatter.hpp
    src/boost/locale/icu/formatters_cache.cpp
    src/boost/locale/icu/formatters_cache.hpp
    src/boost/locale/icu/icu_backend.cpp
    src/boost/locale/icu/icu_backend.hpp
    src/boost/locale/icu/icu_util.hpp
    src/boost/locale/icu/numeric.cpp
    src/boost/locale/icu/time_zone.cpp
    src/boost/locale/icu/time_zone.hpp
    src/boost/locale/icu/uconv.hpp
  )

endif()

if(BOOST_LOCALE_ENABLE_STD)

  target_sources(boost_locale PRIVATE
    src/boost/locale/std/all_generator.hpp
    src/boost/locale/std/codecvt.cpp
    src/boost/locale/std/collate.cpp
    src/boost/locale/std/converter.cpp
    src/boost/locale/std/numeric.cpp
    src/boost/locale/std/std_backend.cpp
    src/boost/locale/std/std_backend.hpp
  )

else()

  target_compile_definitions(boost_locale PRIVATE BOOST_LOCALE_NO_STD_BACKEND=1)

endif()

if(BOOST_LOCALE_ENABLE_WINAPI)

  target_sources(boost_locale PRIVATE
    src/boost/locale/win32/all_generator.hpp
    src/boost/locale/win32/api.hpp
    src/boost/locale/win32/collate.cpp
    src/boost/locale/win32/converter.cpp
    src/boost/locale/win32/numeric.cpp
    src/boost/locale/win32/win_backend.cpp
    src/boost/locale/win32/win_backend.hpp
  )

else()

  target_compile_definitions(boost_locale PRIVATE BOOST_LOCALE_NO_WINAPI_BACKEND=1)

endif()

if(BOOST_LOCALE_ENABLE_WINAPI OR (BOOST_LOCALE_ENABLE_STD AND WIN32))

  target_sources(boost_locale PRIVATE
    src/boost/locale/win32/lcid.cpp
    src/boost/locale/win32/lcid.hpp
  )

endif()

if(BOOST_LOCALE_ENABLE_POSIX)

  target_sources(boost_locale PRIVATE
    src/boost/locale/posix/all_generator.hpp
    src/boost/locale/posix/codecvt.cpp
    src/boost/locale/posix/collate.cpp
    src/boost/locale/posix/converter.cpp
    src/boost/locale/posix/numeric.cpp
    src/boost/locale/posix/posix_backend.cpp
    src/boost/locale/posix/posix_backend.hpp
  )

else()

  target_compile_definitions(boost_locale PRIVATE BOOST_LOCALE_NO_POSIX_BACKEND=1)

endif()

if(BOOST_LOCALE_ENABLE_WINAPI OR BOOST_LOCALE_ENABLE_STD OR BOOST_LOCALE_ENABLE_POSIX)

  target_sources(boost_locale PRIVATE
    src/boost/locale/util/gregorian.cpp
    src/boost/locale/util/gregorian.hpp
  )

endif()

get_target_property(sources boost_locale SOURCES)
set(sources_src ${sources})
list(FILTER sources_src INCLUDE REGEX "^src/")
set(sources_inc ${sources})
list(FILTER sources_inc EXCLUDE REGEX "^src/")
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/src/boost/locale" PREFIX "sources" FILES ${sources_src})
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/include/boost/locale" PREFIX "headers" FILES ${sources_inc})
source_group("headers" FILES "include/boost/locale.hpp")

message(STATUS "Boost.Locale: "
  "iconv ${BOOST_LOCALE_ENABLE_ICONV}, "
  "ICU ${BOOST_LOCALE_ENABLE_ICU}, "
  "POSIX ${BOOST_LOCALE_ENABLE_POSIX}, "
  "std ${BOOST_LOCALE_ENABLE_STD}, "
  "winapi ${BOOST_LOCALE_ENABLE_WINAPI}"
)

option(BOOST_LOCALE_WERROR "Treat warnings as errors in Boost.Locale" ${BOOST_LOCALE_WERROR})
if(MSVC)
    set(BOOST_LOCALE_WARNING_OPTIONS /W4)
    if(BOOST_LOCALE_WERROR)
        list(APPEND BOOST_LOCALE_WARNING_OPTIONS /WX)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(BOOST_LOCALE_WARNING_OPTIONS -Wall -Wextra -Wno-error=deprecated-declarations)
    if(BOOST_LOCALE_WERROR)
        list(APPEND BOOST_LOCALE_WARNING_OPTIONS -Werror)
    endif()
elseif(NOT DEFINED BOOST_LOCALE_WARNING_OPTIONS)
    set(BOOST_LOCALE_WARNING_OPTIONS "")
endif()
target_compile_options(boost_locale PRIVATE ${BOOST_LOCALE_WARNING_OPTIONS})

option(BOOST_LOCALE_ENABLE_CXX11_CHAR_TYPES "Enable C++11 char16/char32 type support (experimental)" OFF)
if(BOOST_LOCALE_ENABLE_CXX11_CHAR_TYPES)
    target_compile_definitions(boost_locale PUBLIC BOOST_LOCALE_ENABLE_CHAR16_T BOOST_LOCALE_ENABLE_CHAR32_T)
endif()

# Testing

if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")

  add_subdirectory(test)
  add_subdirectory(examples)

endif()
