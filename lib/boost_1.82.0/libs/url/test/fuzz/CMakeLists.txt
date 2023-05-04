#
# Copyright (c) 2023 alandefreitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#

source_group("" FILES
        fuzz_parse.cpp
        )

function(add_boost_url_fuzzer NAME SOURCE)
    add_library(fuzzerlib_${NAME} ${SOURCE})
    set_property(TARGET fuzzerlib_${NAME} PROPERTY FOLDER "fuzzing")
    target_link_libraries(fuzzerlib_${NAME} PRIVATE Boost::url)

    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        add_executable(fuzzer_${NAME} ${SOURCE})
        set_property(TARGET fuzzer_${NAME} PROPERTY FOLDER "fuzzing")
        target_link_libraries(fuzzer_${NAME} PRIVATE Boost::url)
        target_compile_options(fuzzer_${NAME} PRIVATE -g -O1 -fsanitize=fuzzer,address,undefined)
        target_link_options(fuzzer_${NAME} PRIVATE -fsanitize=fuzzer,address,undefined)

        add_custom_target(fuzz_${NAME} fuzzer_${NAME} -rss_limit_mb=8192 -max_total_time=30 -timeout=30 DEPENDS fuzz_${NAME})
    endif ()
endfunction()

add_boost_url_fuzzer(parse fuzz_parse.cpp)
