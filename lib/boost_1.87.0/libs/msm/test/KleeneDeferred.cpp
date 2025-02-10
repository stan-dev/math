// Copyright 2010 Christophe Henry
// henry UNDERSCORE christophe AT hotmail DOT com
// This is an extended version of the state machine available in the boost::mpl library
// Distributed under the same license as the original.
// Copyright for the original version:
// Copyright 2005 David Abrahams and Aleksey Gurtovoy. Distributed
// under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/include/std_tuple.hpp>

// back-end
#include <boost/msm/back/state_machine.hpp>
//front-end
#include <boost/msm/front/state_machine_def.hpp>
// functors
#include <boost/msm/front/functor_row.hpp>
// for And_ operator
#include <boost/msm/front/operator.hpp>
#ifndef BOOST_MSM_NONSTANDALONE_TEST
#define BOOST_TEST_MODULE kleene_deferred_test
#endif
#include <boost/test/unit_test.hpp>

using namespace std;
namespace msm = boost::msm;
namespace mpl = boost::mpl;
using namespace msm::front;

namespace
{
    // events
    struct event1 {};
    struct event2 {};


    // front-end: define the FSM structure 
    struct fsm_ : public msm::front::state_machine_def<fsm_>
    {
        // we want deferred events and no state requires deferred events (only the fsm in the
        // transition table), so the fsm does.
        typedef int activate_deferred_events;

        // The list of FSM states
        struct StateA : public msm::front::state<> 
        {
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter=0;
            int exit_counter = 0;
        };
        struct StateB : public msm::front::state<> 
        { 
            template <class Event,class FSM>
            void on_entry(Event const&,FSM& ) {++entry_counter;}
            template <class Event,class FSM>
            void on_exit(Event const&,FSM& ) {++exit_counter;}
            int entry_counter = 0;
            int exit_counter = 0;
        };


        // the initial state of the player SM. Must be defined
        typedef StateA initial_state;

        struct is_event1
        {
            template <class EVT, class FSM, class SourceState, class TargetState>
            bool operator()(EVT const& evt, FSM&, SourceState&, TargetState&)
            {
                bool is_deferred = boost::any_cast<event1>(&evt) != 0;
                return is_deferred;
            }
        };

        typedef fsm_ p; // makes transition table cleaner

        // Transition table for player
        struct transition_table : boost::fusion::vector<
            //    Start     Event         Next      Action               Guard
            //  +---------+-------------+---------+---------------------+----------------------+
            Row < StateA  , event1      , StateB  , none                , none                 >,
            Row < StateB  , boost::any  , none    , Defer               , is_event1            >,
            Row < StateB  , event2      , StateA  , none                , none           >
            //  +---------+-------------+---------+---------------------+----------------------+
        >/*>::type*/ {};
        // Replaces the default no-transition response.
        template <class FSM,class Event>
        void no_transition(Event const&, FSM&,int)
        {
            BOOST_FAIL("no_transition called!");
        }
        
    };
    // Pick a back-end
    typedef msm::back::state_machine<fsm_> Fsm;

    BOOST_AUTO_TEST_CASE(kleene_deferred_test)
    {     
        Fsm fsm;

        fsm.start(); 
        BOOST_CHECK_MESSAGE(fsm.get_state<fsm_::StateA&>().entry_counter == 1,"StateA entry not called correctly");

        fsm.process_event(event1()); 
        BOOST_CHECK_MESSAGE(fsm.current_state()[0] == 1,"StateB should be active"); 
        BOOST_CHECK_MESSAGE(fsm.get_state<fsm_::StateA&>().exit_counter == 1,"StateA exit not called correctly");
        BOOST_CHECK_MESSAGE(fsm.get_state<fsm_::StateB&>().entry_counter == 1,"StateB entry not called correctly");

        fsm.process_event(event1());
        BOOST_CHECK_MESSAGE(fsm.current_state()[0] == 1, "StateB should be active");
        BOOST_CHECK_MESSAGE(fsm.get_state<fsm_::StateA&>().exit_counter == 1, "StateA exit not called correctly");
        BOOST_CHECK_MESSAGE(fsm.get_state<fsm_::StateB&>().entry_counter == 1, "StateB entry not called correctly");

        fsm.process_event(event2());
        BOOST_CHECK_MESSAGE(fsm.current_state()[0] == 1, "StateB should be active");
        BOOST_CHECK_MESSAGE(fsm.get_state<fsm_::StateA&>().exit_counter == 2, "StateA exit not called correctly");
        BOOST_CHECK_MESSAGE(fsm.get_state<fsm_::StateB&>().entry_counter == 2, "StateB entry not called correctly");

    }
}

