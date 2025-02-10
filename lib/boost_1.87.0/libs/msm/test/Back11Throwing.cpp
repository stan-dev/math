// Copyright 2024 Christophe Henry
// henry UNDERSCORE christophe AT hotmail DOT com
// This is an extended version of the state machine available in the boost::mpl library
// Distributed under the same license as the original.
// Copyright for the original version:
// Copyright 2005 David Abrahams and Aleksey Gurtovoy. Distributed
// under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// back-end
#include <boost/msm/back11/state_machine.hpp>
//front-end
#include <boost/msm/front/state_machine_def.hpp>
#include <boost/msm/front/functor_row.hpp>
#include <string>
#include <iostream>

#ifndef BOOST_MSM_NONSTANDALONE_TEST
#define BOOST_TEST_MODULE Back11ThrowingTest
#endif
#include <boost/test/unit_test.hpp>

namespace msm = boost::msm;
namespace mpl = boost::mpl;
using namespace boost::msm::front;

namespace
{
    template <typename T>
    void print(T t)
    {
        std::cout << t << "\n";
    }

    template <typename... Ts>
    using Row = boost::msm::front::Row<Ts...>;

    using None = boost::msm::front::none;

    // events
    struct ExceptionRaised {};

    // states
    struct Init : boost::msm::front::state<> 
    {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter;}
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter;}
        int entry_counter = 0;
        int exit_counter = 0;
    };
    struct WaitingForException : boost::msm::front::state<> 
    {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter;}
        int entry_counter = 0;
        int exit_counter = 0;
    };
    struct Intermediate : boost::msm::front::state<> 
    {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter;}
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter;}
        int entry_counter = 0;
        int exit_counter = 0;
    };
    struct End : boost::msm::front::terminate_state<> 
    {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter;}
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter;}
        int entry_counter = 0;
        int exit_counter = 0;
    };

    // actions
    struct ThrowingAction
    {
        template <class Event, class StateMachine, class SourceState, class TargetState>
        void operator()(const Event&, StateMachine&, SourceState&, TargetState&)
        {
            throw std::runtime_error("foo");
        }
    };
    struct ExceptionHandler
    {
        template <class Event, class StateMachine, class SourceState, class TargetState>
        void operator()(const Event&, StateMachine&, SourceState&, TargetState&)
        {
        }
    };

    struct MyMachineFrontend : public boost::msm::front::state_machine_def<MyMachineFrontend>
    {
        using initial_state = boost::mpl::vector<Init,WaitingForException>;

        using transition_table = boost::mpl::vector<
            Row<Init, None, Intermediate, ThrowingAction>,

            Row<WaitingForException, ExceptionRaised, End, ExceptionHandler>
        >;

        template <class StateMachine, class Event>
        void exception_caught(const Event&, StateMachine& machine, std::exception&)
        {
            machine.process_event(ExceptionRaised{});
        }
    };

    using MyMachine = boost::msm::back11::state_machine<MyMachineFrontend>;


    BOOST_AUTO_TEST_CASE(Back11ThrowingTest)
    {
        MyMachine m;
        m.start();

        BOOST_CHECK_MESSAGE(m.get_state<Init&>().entry_counter == 1, "Init entry not called correctly");
        BOOST_CHECK_MESSAGE(m.get_state<Init&>().exit_counter == 1, "Init exit not called correctly");
        BOOST_CHECK_MESSAGE(m.get_state<Intermediate&>().entry_counter == 0, "Intermediate entry not called correctly");
        BOOST_CHECK_MESSAGE(m.get_state<Intermediate&>().exit_counter == 0, "Intermediate exit not called correctly");

        BOOST_CHECK_MESSAGE(m.get_state<WaitingForException&>().entry_counter == 1, "WaitingForException entry not called correctly");
        BOOST_CHECK_MESSAGE(m.get_state<WaitingForException&>().exit_counter == 1, "WaitingForException exit not called correctly");
        BOOST_CHECK_MESSAGE(m.get_state<End&>().entry_counter == 1, "End entry not called correctly");
        BOOST_CHECK_MESSAGE(m.get_state<End&>().exit_counter == 0, "End exit not called correctly");


        BOOST_CHECK_MESSAGE(m.current_state()[0] == 0, "Init should be active");
        BOOST_CHECK_MESSAGE(m.current_state()[1] == 3, "End should be active");
    }
}


