if(NOT TARGET tests)
    add_custom_target(tests)
    set_property(TARGET tests PROPERTY FOLDER _deps)
endif()

add_library(boost_cobalt_static_tests concepts.cpp util.cpp)
target_link_libraries(boost_cobalt_static_tests Boost::cobalt)

add_executable(boost_cobalt_main         main.cpp)
add_executable(boost_cobalt_main_compile main_compile.cpp)
add_executable(boost_cobalt_basic_tests
      async_for.cpp test_main.cpp promise.cpp with.cpp op.cpp handler.cpp join.cpp race.cpp this_coro.cpp leaf.cpp
      channel.cpp generator.cpp run.cpp task.cpp gather.cpp wait_group.cpp wrappers.cpp left_race.cpp
      strand.cpp fork.cpp thread.cpp any_completion_handler.cpp detached.cpp monotonic_resource.cpp sbo_resource.cpp)

target_link_libraries(boost_cobalt_main         Boost::cobalt)
target_link_libraries(boost_cobalt_main_compile Boost::cobalt)
target_link_libraries(boost_cobalt_basic_tests  Boost::cobalt Boost::unit_test_framework)

add_test(NAME boost_cobalt_main COMMAND boost_cobalt_main)
add_test(NAME boost_cobalt_basic_tests COMMAND boost_cobalt_basic_tests)
target_compile_features(boost_cobalt PUBLIC cxx_std_20)

add_dependencies(tests boost_cobalt_main boost_cobalt_basic_tests boost_cobalt_static_tests)
