enable_testing()

add_library(boost_process_v2_test_impl test_impl.cpp)
target_compile_definitions(boost_process_v2_test_impl PUBLIC -DBOOST_PROCESS_V2_SEPARATE_COMPILATION=1)

if (WIN32)
    target_compile_definitions(boost_process_v2_test_impl PUBLIC WIN32_LEAN_AND_MEAN=1)
    target_link_libraries(boost_process_v2_test_impl Boost::process Boost::unit_test_framework Boost::process Ntdll)
else()
    target_link_libraries(boost_process_v2_test_impl Boost::process Boost::unit_test_framework Boost::process)
endif()

function(boost_process_v2_standalone_test name)
    add_executable(boost_process_v2_${name} ${name}.cpp)
    target_link_libraries(boost_process_v2_${name} Boost::process Boost::system Boost::filesystem boost_process_v2_test_impl)
    add_test(NAME boost_process_v2_${name} COMMAND $<TARGET_FILE:boost_process_v2_${name}> )
endfunction()

boost_process_v2_standalone_test(utf8)
boost_process_v2_standalone_test(cstring_ref)
boost_process_v2_standalone_test(pid)
boost_process_v2_standalone_test(environment)
boost_process_v2_standalone_test(shell)

add_executable(boost_process_v2_test_target target.cpp)
target_link_libraries(boost_process_v2_test_target PUBLIC Boost::process Boost::system)

function(boost_process_v2_test_with_target name)
    add_executable(boost_process_v2_${name} ${name}.cpp)
    target_link_libraries(boost_process_v2_${name} Boost::process Boost::system Boost::filesystem boost_process_v2_test_impl)
    add_dependencies(boost_process_v2_${name} boost_process_v2_test_target)
    add_test(NAME boost_process_v2_${name} COMMAND $<TARGET_FILE:boost_process_v2_${name}>
                                                -- $<TARGET_FILE:boost_process_v2_test_target>)


endfunction()

boost_process_v2_test_with_target(process)
boost_process_v2_test_with_target(ext)