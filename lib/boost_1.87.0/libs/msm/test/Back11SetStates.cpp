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
#include <boost/msm/back11/state_machine.hpp>
//front-end
#include <boost/msm/front/state_machine_def.hpp>
#include <boost/msm/front/functor_row.hpp>

#ifndef BOOST_MSM_NONSTANDALONE_TEST
#define BOOST_TEST_MODULE back11_set_states_test
#endif
#include <boost/test/unit_test.hpp>

namespace msm = boost::msm;
namespace mpl = boost::mpl;
namespace fusion = boost::fusion;
using namespace msm::front;

namespace
{
    struct InputEvent {};
    struct ExitEvent
    {
        ExitEvent() = default;
        template <typename T>
        ExitEvent(const T&) {};
    };

    struct State : public boost::msm::front::state<>
    {
        template <class Event, class FSM>
        void on_entry(Event const& e, FSM& fsm)
        {
            std::cout << "State - on_entry" << "\n";
        };

        template <class Event, class FSM>
        void on_exit(Event const& e, FSM& fsm)
        {
            std::cout << "State - on_exit" << "\n";
        }
    };

    struct TerminateState : public boost::msm::front::terminate_state<>
    {
        template <class Event, class FSM>
        void on_entry(Event const& e, FSM& fsm)
        {
            std::cout << "TerminateState - on_entry" << "\n";
        };

        template <class Event, class FSM>
        void on_exit(Event const& e, FSM& fsm)
        {
            std::cout << "TerminateState - on_exit" << "\n";
        }
    };

    struct PseudoExitState : public boost::msm::front::exit_pseudo_state<ExitEvent>
    {
    };

    struct SubFsm_ : public msm::front::state_machine_def<SubFsm_>
    {
        template <class Event, class FSM>
        void on_entry(Event const& e, FSM& fsm)
        {
            std::cout << "SubFsm - on_entry" << "\n";
        };

        template <class Event, class FSM>
        void on_exit(Event const& e, FSM& fsm)
        {
            std::cout << "SubFsm - on_exit" << "\n";
        }

        typedef State initial_state;

        struct transition_table : mpl::vector<
            Row< State, InputEvent, PseudoExitState, none, none >
        > {};

        template <class FSM, class Event>
        void no_transition(Event const& e, FSM&, int state)
        {
            std::cout << "SubFsm: no transition from state " << state
                << " on event " << typeid(e).name() << std::endl;
        }
    };
    typedef msm::back11::state_machine<SubFsm_> SubFsm;

    struct Fsm_ : public msm::front::state_machine_def<Fsm_>
    {
        template <class Event, class FSM>
        void on_entry(Event const& e, FSM& fsm)
        {
            std::cout << "Fsm - on_entry" << "\n";
        };

        template <class Event, class FSM>
        void on_exit(Event const& e, FSM& fsm)
        {
            std::cout << "Fsm - on_exit" << "\n";
        }

        typedef SubFsm initial_state;

        struct transition_table : mpl::vector<
            Row< SubFsm::exit_pt<PseudoExitState>, ExitEvent, TerminateState, none, none >
        > {};

        template <class FSM, class Event>
        void no_transition(Event const& e, FSM&, int state)
        {
            std::cout << "Fsm: no transition from state " << state
                << " on event " << typeid(e).name() << std::endl;
        }
    };
    typedef msm::back11::state_machine<Fsm_> Fsm;


    BOOST_AUTO_TEST_CASE(back11_set_states_test)
    {     
        InputEvent input_event;

        std::cout << "\nSubFsms created in constructor\n";
        Fsm fsm_subfsms_created_in_constructor(boost::msm::back::states_
            << SubFsm());
        fsm_subfsms_created_in_constructor.start();
        fsm_subfsms_created_in_constructor.process_event(input_event);

        std::cout << "\nSubFsms created through set_states method\n";
        Fsm fsm_subfsms_created_thorugh_constructor;
        fsm_subfsms_created_thorugh_constructor.set_states(boost::msm::back::states_
            << SubFsm());
        fsm_subfsms_created_thorugh_constructor.start();
        fsm_subfsms_created_thorugh_constructor.process_event(input_event);
    }
}

