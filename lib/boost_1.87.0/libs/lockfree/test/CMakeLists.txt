if(NOT TARGET boost_lockfree_all_tests)
    add_custom_target(boost_lockfree_all_tests)
endif()

include(BoostTest)

if (BOOST_LOCKFREE_BUILD_TESTS)
    set(BUILD_TESTING TRUE)
endif()

if (NOT BUILD_TESTING AND NOT TARGET tests)
    add_custom_target(tests)
endif()

add_library(boost_lockfree_test_common INTERFACE)

if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.23)
  target_sources(boost_lockfree_test_common PUBLIC FILE_SET HEADERS FILES test_common.hpp test_helpers.hpp )
  source_group( TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES test_common.hpp test_helpers.hpp)
endif()

set(Tests
    destructor_test
    freelist_test
    queue_bounded_stress_test
    queue_fixedsize_stress_test
    queue_interprocess_test
    queue_test
    queue_unbounded_stress_test
    spsc_queue_stress_test
    spsc_queue_test
    stack_bounded_stress_test
    stack_fixedsize_stress_test
    stack_interprocess_test
    stack_test
    stack_unbounded_stress_test
    tagged_ptr_test
)

foreach(Test ${Tests})
    set (Libs Boost::lockfree
              Boost::unit_test_framework
              Boost::thread
              boost_lockfree_test_common)

    if (Test MATCHES ".*interprocess.*")
        list(APPEND Libs Boost::interprocess)
    endif()

    source_group( TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${Test}.cpp)

    boost_test(
        TYPE run
        PREFIX boost_lockfree
        NAME ${Test}
        SOURCES ${Test}.cpp
        LINK_LIBRARIES ${Libs}
        COMPILE_DEFINITIONS BOOST_TEST_NO_OLD_TOOLS
    )
    set_target_properties( boost_lockfree-${Test} PROPERTIES CMAKE_CXX_STANDARD_REQUIRED 14)

    if (TARGET boost_lockfree-${Test})
        add_dependencies(boost_lockfree_all_tests boost_lockfree-${Test} )
    endif()

endforeach()
