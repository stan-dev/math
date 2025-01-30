// Copyright 2024 Christophe Henry
// henry UNDERSCORE christophe AT hotmail DOT com
// This is an extended version of the state machine available in the boost::mpl library
// Distributed under the same license as the original.
// Copyright for the original version:
// Copyright 2005 David Abrahams and Aleksey Gurtovoy. Distributed
// under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>

#ifdef __APPLE__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#endif
#include <boost/msm/back/state_machine.hpp>
#include <boost/msm/front/euml/common.hpp>
#include <boost/msm/front/functor_row.hpp>
#include <boost/msm/front/state_machine_def.hpp>
#ifdef __APPLE__
#pragma clang diagnostic pop
#endif

#ifndef BOOST_MSM_NONSTANDALONE_TEST
#define BOOST_TEST_MODULE back11_anonymous_and_guard_test
#endif
#include <boost/test/unit_test.hpp>


namespace bmf = boost::msm::front;

// events
struct success {
    success(size_t value) : value(value) {}
    size_t value;
};
struct failure {
};

// frontend
struct Bug : public bmf::state_machine_def<Bug> {
    // list of FSM states
    struct Init : public bmf::state<> {
    };
    struct Idle : public bmf::state<> {
    };
    struct Completed : public bmf::terminate_state<> {
    };

    // defines an orthogonal region
    struct Running : public bmf::state<> {
    };

    // actions
    struct log_action {
        template<class Event, class FSM, class SourceState, class TargetState>
        void operator()(const Event&, FSM&, SourceState&, TargetState&)
        {
        }
    };

    // guards
    struct validate {
        template<class Event, class FSM, class SourceState, class TargetState>
        bool operator()(const Event& e, FSM&, SourceState&, TargetState&)
        {
            return e.value > 0;
        }
    };

    // state machine properties
    typedef boost::mpl::vector2<Running, Init> initial_state;

    // Transition table for traceroute
    struct transition_table : boost::mpl::vector<
        //    Start      Event           Next      Action                     Guard
        //  +----------+---------------+-----------+------------------------+--------------------+
        bmf::Row< Running, success, bmf::none, log_action, validate           >,
        //  +---------+-------------+---------+---------------------------+----------------------+
        bmf::Row< Init, success, Idle, log_action, bmf::none          >,
        bmf::Row< Init, failure, Idle, log_action, bmf::none          >,
        bmf::Row< Idle, bmf::none, Completed, log_action, bmf::none          >
        //  +---------+-------------+---------+---------------------------+----------------------+
    > {};
};

// backend
typedef boost::msm::back::state_machine<Bug> MyStateMachine;

BOOST_AUTO_TEST_CASE(back11_anonymous_and_guard_test1)
{
    MyStateMachine sm;
    sm.start();
    sm.process_event(success(0));
    BOOST_CHECK_MESSAGE(sm.current_state()[0] == 0, "Running should be active");
    BOOST_CHECK_MESSAGE(sm.current_state()[1] == 3, "Completed should be active");
}

BOOST_AUTO_TEST_CASE(back11_anonymous_and_guard_test2)
{
    MyStateMachine sm;
    sm.start();
    sm.process_event(success(1));
    BOOST_CHECK_MESSAGE(sm.current_state()[0] == 0, "Running should be active");
    BOOST_CHECK_MESSAGE(sm.current_state()[1] == 3, "Completed should be active");
}
