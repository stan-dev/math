// Copyright 2010 Christophe Henry
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
#include <boost/msm/front/state_machine_def.hpp>

#include <boost/test/unit_test.hpp>

namespace msm = boost::msm;
namespace mpl = boost::mpl;
using namespace boost::msm::front;

namespace
{
    // events
    struct event1 {};
    struct event2 {};


    // front-end: define the FSM structure
    struct my_machine_ : public msm::front::state_machine_def<my_machine_>
    {

        my_machine_()
        {}

        // The list of FSM states
        struct State1 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State2 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State3 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State4 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State5 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State6 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State7 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State8 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State9 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State10 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State11 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };
        struct State12 : public msm::front::state<>
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter=0;
        };

        // the initial state of the player SM. Must be defined
        typedef State1 initial_state;

        // transition actions

        // guard conditions

        typedef my_machine_ p; // makes transition table cleaner

        // Transition table for player
        struct transition_table : mpl::vector<
            //    Start     Event         Next      Action               Guard
            //  +---------+-------------+---------+---------------------+----------------------+
            _row < State1  , event1      , State2                                               >,
            _row < State2  , event1      , State3                                               >,
            _row < State3  , event1      , State4                                               >,
            _row < State4  , event1      , State5                                               >,
            _row < State5  , event1      , State6                                               >,
            _row < State6  , event1      , State7                                               >,
            _row < State7  , event1      , State8                                               >,
            _row < State8  , event1      , State9                                               >,
            _row < State9  , event1      , State10                                              >,
            _row < State10 , event1      , State11                                              >,
            _row < State11 , event1      , State12                                              >,

            _row < State12 , event2      , State11                                              >,
            _row < State11 , event2      , State10                                              >,
            _row < State10 , event2      , State9                                               >,
            _row < State9  , event2      , State8                                               >,
            _row < State8  , event2      , State7                                               >,
            _row < State7  , event2      , State6                                               >,
            _row < State6  , event2      , State5                                               >,
            _row < State5  , event2      , State4                                               >,
            _row < State4  , event2      , State3                                               >,
            _row < State3  , event2      , State2                                               >,
            _row < State2  , event2      , State1                                               >
            //  +---------+-------------+---------+---------------------+----------------------+
        > {};
        // Replaces the default no-transition response.
        template <class FSM,class Event>
        void no_transition(Event const&, FSM&,int)
        {
            BOOST_FAIL("no_transition called!");
        }
        // init counters
        template <class Event,class FSM>
        void on_entry(Event const&,FSM&)
        {

        }
    };
    // Pick a back-end
    typedef msm::back::state_machine<my_machine_> my_machine;


    BOOST_AUTO_TEST_CASE( big_with_functors_test )
    {
        my_machine p;

        // needed to start the highest-level SM.
        p.start();
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 0,"State1 should be active"); //State1
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State1&>().exit_counter == 0,"State1 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State1&>().entry_counter == 1,"State1 entry not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State5&>().exit_counter == 0,"State5 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State5&>().entry_counter == 0,"State5 entry not called correctly");

        p.process_event(event1());
        p.process_event(event1());
        p.process_event(event1());
        p.process_event(event1());
        p.process_event(event1());
        p.process_event(event1());
        p.process_event(event1());
        p.process_event(event1());
        p.process_event(event1());
        p.process_event(event1());
        p.process_event(event1());

        BOOST_CHECK_MESSAGE(p.current_state()[0] == 11,"State12 should be active"); //State12
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State12&>().exit_counter == 0,"State12 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State12&>().entry_counter == 1,"State12 entry not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State2&>().exit_counter == 1,"State2 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State2&>().entry_counter == 1,"State2 entry not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State7&>().exit_counter == 1,"State7 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State7&>().entry_counter == 1,"State7 entry not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State10&>().exit_counter == 1,"State10 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State10&>().entry_counter == 1,"State10 entry not called correctly");


        p.process_event(event2());
        p.process_event(event2());
        p.process_event(event2());
        p.process_event(event2());
        p.process_event(event2());
        p.process_event(event2());
        p.process_event(event2());
        p.process_event(event2());
        p.process_event(event2());
        p.process_event(event2());
        p.process_event(event2());

        BOOST_CHECK_MESSAGE(p.current_state()[0] == 0,"State1 should be active"); //State1
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State1&>().exit_counter == 1,"State1 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State1&>().entry_counter == 2,"State1 entry not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State2&>().exit_counter == 2,"State2 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State2&>().entry_counter == 2,"State2 entry not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State7&>().exit_counter == 2,"State7 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State7&>().entry_counter == 2,"State7 entry not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State10&>().exit_counter == 2,"State10 exit not called correctly");
        BOOST_CHECK_MESSAGE(p.get_state<my_machine_::State10&>().entry_counter == 2,"State10 entry not called correctly");

    }
}

