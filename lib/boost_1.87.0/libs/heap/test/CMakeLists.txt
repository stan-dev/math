if(NOT TARGET boost_heap_all_tests)
    add_custom_target(boost_heap_all_tests)
endif()

add_library(boost_heap_test_common INTERFACE)

include(BoostTest)

if (BOOST_HEAP_BUILD_TESTS)
    set(BUILD_TESTING TRUE)
endif()


if (NOT BUILD_TESTING AND NOT TARGET tests)
    add_custom_target(tests)
endif()


if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.23)
  target_sources(boost_heap_test_common PUBLIC FILE_SET HEADERS FILES
      common_heap_tests.hpp
      merge_heap_tests.hpp
      mutable_heap_tests.hpp
      stable_heap_tests.hpp
  )
  source_group( TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES
      common_heap_tests.hpp
      merge_heap_tests.hpp
      mutable_heap_tests.hpp
      stable_heap_tests.hpp)
endif()

set(Tests
    binomial_heap_test
    d_ary_heap_test
    fibonacci_heap_test
    pairing_heap_tests
    priority_queue_test
    skew_heap_test
    mutable_heap_test
)

foreach(Test ${Tests})
    set (Libs Boost::heap
              Boost::unit_test_framework
              Boost::container
              Boost::foreach
              boost_heap_test_common)

    source_group( TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${Test}.cpp)

    boost_test(
        TYPE run
        PREFIX boost_heap
        NAME ${Test}
        SOURCES ${Test}.cpp
        LINK_LIBRARIES ${Libs}
    )

    if (TARGET boost_heap-${Test})
        add_dependencies(boost_heap_all_tests boost_heap-${Test} )
    endif()
endforeach()
