// Copyright 2010 Christophe Henry
// henry UNDERSCORE christophe AT hotmail DOT com
// This is an extended version of the state machine available in the boost::mpl library
// Distributed under the same license as the original.
// Copyright for the original version:
// Copyright 2005 David Abrahams and Aleksey Gurtovoy. Distributed
// under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>

#ifndef BOOST_MSM_NONSTANDALONE_TEST
#define BOOST_TEST_MODULE back11_transition_skipping
#endif
#include <boost/test/unit_test.hpp>
#include <boost/core/demangle.hpp>
//#define USE_PRE11_MSM_BACKEND   // <- toggle this to switch between back vs back11

// back-end
#if defined(USE_PRE11_MSM_BACKEND)
#include <boost/msm/back/state_machine.hpp>
namespace msm_back = boost::msm::back;
#else
#include <boost/msm/back11/state_machine.hpp>
namespace msm_back = boost::msm::back11;
#endif

// front-end
#include <boost/fusion/include/insert_range.hpp>
#include <boost/msm/front/euml/operator.hpp>
#include <boost/msm/front/functor_row.hpp>
#include <boost/msm/front/state_machine_def.hpp>


template <typename ...T>
using msm_flags = boost::mpl::vector<T...>;
template <typename ...T>
using msm_transition_table = boost::mpl::vector<T...>;



namespace msmf = boost::msm::front;
namespace euml = boost::msm::front::euml;

struct evt1 {};
struct evt2 {};
struct evt3 {};

struct top_ : public msmf::state_machine_def<top_> {
    int count_guard_false = 0;
    int count_guard_true = 0;


    struct nested : public msmf::state<> {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter=0;
        int exit_counter=0;
    };

    struct other1 : public msmf::state<> {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
    };

    struct other2 : public msmf::state<> {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
    };

    struct other3 : public msmf::state<> {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
    };

    struct other4 : public msmf::state<> {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
    };
    struct other5 : public msmf::state<> {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
    };

    struct guard_true {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM& fsm, SourceState&, TargetState&) const {
            ++fsm.count_guard_true;
            return true;
        }
    };

    struct guard_false {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM& fsm, SourceState&, TargetState&) const {
            ++fsm.count_guard_false;
            return false;
        }
    };


    typedef nested initial_state;

    using transition_table = msm_transition_table<
        //          Start        Event     Next         Action      Guard
        //--------+------------+---------+------------+-----------+--------------
        msmf::Row < nested, msmf::none, other1, msmf::none, guard_false>,
        msmf::Row < nested, msmf::none, other2, msmf::none, guard_true >,  // <- this transition is not executed with back11
        msmf::Row < nested, msmf::none, other3, msmf::none, guard_false>,
        msmf::Row < nested, msmf::none, other4, msmf::none, guard_false>,
        msmf::Row < other1, evt3, nested, msmf::none, msmf::none       >
    >;
};

using top = msm_back::state_machine<top_>;
   

BOOST_AUTO_TEST_CASE(back11_transition_skipping)
{
    top msm;
    msm.start();
    BOOST_CHECK_MESSAGE(msm.count_guard_true == 1, "guard_true should be called once");
    BOOST_CHECK_MESSAGE(msm.count_guard_false == 2, "guard_false should be called once");

    BOOST_CHECK_MESSAGE(msm.get_state<top_::nested&>().exit_counter == 1, "nested exit not called correctly");
    BOOST_CHECK_MESSAGE(msm.get_state<top_::nested&>().entry_counter == 1, "nested entry not called correctly");
    BOOST_CHECK_MESSAGE(msm.get_state<top_::other2&>().exit_counter == 0, "other2 exit not called correctly");
    BOOST_CHECK_MESSAGE(msm.get_state<top_::other2&>().entry_counter == 1, "other2 entry not called correctly");
    BOOST_CHECK_MESSAGE(msm.get_state<top_::other1&>().exit_counter == 0, "other1 exit not called correctly");
    BOOST_CHECK_MESSAGE(msm.get_state<top_::other1&>().entry_counter == 0, "other1 entry not called correctly");
    BOOST_CHECK_MESSAGE(msm.get_state<top_::other3&>().exit_counter == 0, "other3 exit not called correctly");
    BOOST_CHECK_MESSAGE(msm.get_state<top_::other3&>().entry_counter == 0, "other3 entry not called correctly");
    BOOST_CHECK_MESSAGE(msm.get_state<top_::other4&>().exit_counter == 0, "other4 exit not called correctly");
    BOOST_CHECK_MESSAGE(msm.get_state<top_::other4&>().entry_counter == 0, "other4 entry not called correctly");
    BOOST_CHECK_MESSAGE(msm.current_state()[0] == 2, "other2 should be active"); //Open

}
