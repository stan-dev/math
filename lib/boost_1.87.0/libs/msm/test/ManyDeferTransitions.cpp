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
#include <boost/msm/back/state_machine.hpp>
//front-end
#include <boost/msm/front/functor_row.hpp>
#include <boost/msm/front/state_machine_def.hpp>
#ifndef BOOST_MSM_NONSTANDALONE_TEST
#define BOOST_TEST_MODULE back_many_deferred_transitions
#endif
#include <boost/test/unit_test.hpp>

namespace msm = boost::msm;
namespace msmf = boost::msm::front;
namespace mpl = boost::mpl;

namespace
{
    // ----- Events
    struct Event1 {};
    struct Event2 {};
    struct Event3 {};
    struct Event4 {};

    // ----- State machine
    struct Sm1_ :msmf::state_machine_def<Sm1_> {

        // States
        struct State1 :msmf::state<> {
            template <class Event, class Fsm>
            void on_entry(Event const&, Fsm&) const {
            }
        };
        struct State2 :msmf::state<> {
            template <class Event, class Fsm>
            void on_entry(Event const&, Fsm&) const {
            }
        };
        struct State3 :msmf::state<> {
            template <class Event, class Fsm>
            void on_entry(Event const&, Fsm&) const {
            }
        };
        struct State4 :msmf::state<> {
            template <class Event, class Fsm>
            void on_entry(Event const&, Fsm&) const {
            }
        };
        struct State5 :msmf::state<> {
            template <class Event, class Fsm>
            void on_entry(Event const&, Fsm&) const {
            }
        };
        // Set initial state
        typedef State1 initial_state;
        // Enable deferred capability
        typedef int activate_deferred_events;
        // Transition table
        struct transition_table :mpl::vector<
            //          Start   Event   Next        Action       Guard
            msmf::Row < State1, Event1, msmf::none, msmf::Defer>,
            msmf::Row < State1, Event2, msmf::none, msmf::Defer>,
            msmf::Row < State1, Event3, msmf::none, msmf::Defer>,
            msmf::Row < State1, Event4, State2    , msmf::none>,

            msmf::Row < State2, Event3, msmf::none, msmf::Defer>,
            msmf::Row < State2, Event1, msmf::none, msmf::Defer>,
            msmf::Row < State2, Event2, State3    , msmf::none>,

            msmf::Row < State3, Event1, State4    , msmf::none>,
            msmf::Row < State3, Event3, State5    , msmf::none>,
            msmf::Row < State4, Event3, State5    , msmf::none>,
            msmf::Row < State5, Event1, State4    , msmf::none>
        > {};

        template <class Fsm, class Event>
        void no_transition(Event const&, Fsm&, int /*state*/) {
        }
    };

    // Pick a back-end
    typedef msm::back::state_machine<Sm1_> Sm1;


    BOOST_AUTO_TEST_CASE(back_many_deferred_transitions)
    {
        Sm1 sm1;
        sm1.start();
        sm1.process_event(Event1());
        sm1.process_event(Event2());
        sm1.process_event(Event3());
        sm1.process_event(Event4());

        BOOST_CHECK_MESSAGE(sm1.current_state()[0] == 4, "State5 should be active");
    }
}

