
file(GLOB_RECURSE ALL_EXAMPLES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

find_package(OpenSSL)

foreach(SRC ${ALL_EXAMPLES})
    get_filename_component(NAME ${SRC} NAME_WLE )
    # ticker requires
    if (NAME STREQUAL ticker)
        if (TARGET Boost::json)
            add_executable(boost_cobalt_example_${NAME} ${SRC} )
            target_link_libraries(boost_cobalt_example_${NAME} PUBLIC Boost::cobalt Boost::json Boost::url OpenSSL::SSL)
            target_compile_definitions(boost_cobalt_example_${NAME} PUBLIC)
        endif()
        continue()
    endif()

    if (NAME STREQUAL python)
        find_package(Python 3.8 COMPONENTS Interpreter Development.Module)
        if (NOT Python_FOUND)
            message(WARNING "Python not found, skipping python example")
            continue()
        endif()
        execute_process(
                COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
        list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
        find_package(nanobind CONFIG)
        if (NOT nanobind_FOUND)
            message(WARNING "nanobind not found, skipping python example")
            continue()
        endif()
        nanobind_add_module(boost_cobalt_example_python python.cpp)
        target_link_libraries(boost_cobalt_example_python PRIVATE Boost::cobalt)
        add_custom_command(
                TARGET boost_cobalt_example_python
                COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/python.py
                ${CMAKE_CURRENT_BINARY_DIR}/python.py)
        continue()
    endif()

    add_executable(boost_cobalt_example_${NAME} ${SRC})
    target_link_libraries(boost_cobalt_example_${NAME} PUBLIC Boost::cobalt OpenSSL::SSL)
    target_compile_definitions(boost_cobalt_example_${NAME} PUBLIC)
endforeach()
