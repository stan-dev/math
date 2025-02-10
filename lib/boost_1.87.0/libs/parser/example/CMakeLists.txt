# Copyright (C) 2020 T. Zachary Laine
#
# Distributed under the Boost Software License, Version 1.0. (See
# accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
include_directories(${CMAKE_HOME_DIRECTORY})

include(CTest)

enable_testing()

add_custom_target(run_examples COMMAND ${CMAKE_CTEST_COMMAND} -VV -C ${CMAKE_CFG_INTDIR})

macro(add_example name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} parser ${ARGN})
    set_property(TARGET ${name} PROPERTY CXX_STANDARD ${CXX_STD})
    add_test(${name} ${CMAKE_CURRENT_BINARY_DIR}/${name})
endmacro()

set(libs)
if (BUILD_WITH_HANA)
  set(libs boost)
endif()

add_example(hello ${libs})
add_example(trivial ${libs})
add_example(trivial_skipper ${libs})
add_example(semantic_actions ${libs})
add_example(rule_intro ${libs})
add_example(struct_rule ${libs})
add_example(parsing_into_a_struct ${libs})
add_example(parsing_into_a_class ${libs})
add_example(user_error_handler ${libs})

if (BUILD_WITH_HANA)
  add_example(roman_numerals boost)
  add_example(self_filling_symbol_table boost)
  add_example(json boost)
  add_example(callback_json boost)
endif()
