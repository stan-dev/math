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
#include <boost/msm/front/puml/puml.hpp>
#include <PumlCommon.hpp>

#ifndef BOOST_MSM_NONSTANDALONE_TEST
#define BOOST_TEST_MODULE back11_only_string_puml_test
#endif
#include <boost/test/unit_test.hpp>

using namespace std;
namespace msm = boost::msm;
using namespace msm::front;
using namespace msm::front::puml;

namespace
{
    // note in the puml tests will be non-specialized types marked with _ 

    // front-end: define the FSM structure 
    struct front_ : public msm::front::state_machine_def<front_>
    {
        BOOST_MSM_PUML_DECLARE_TABLE(
            R"( 
            @startuml Player
            skinparam linetype polyline
            state Player{
                [*]-> StateA
                StateA          -> StateB        : some_event
                StateB          -> StateA        : some_event                
                --
                [*]-> StateC
                StateC          -> TerminalState : terminate_event
                TerminalState   -> [*]
            }
            @enduml
        )"
        )

        // Replaces the default no-transition response.
        template <class FSM,class Event>
        void no_transition(Event const&, FSM&,int)
        {
            BOOST_FAIL("no_transition called!");
        }
    };
    // Pick a back-end
    typedef msm::back11::state_machine<front_> machine;


    BOOST_AUTO_TEST_CASE(back11_string_terminate_puml_test)
    {     
        machine p;
        static_assert(msm::back11::get_number_of_regions<typename machine::initial_state>::type::value == 2);
        static_assert(::boost::mpl::size<typename machine::transition_table>::type::value == 3);

        p.start(); 
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 0, "StateA should be active");
        BOOST_CHECK_MESSAGE(p.current_state()[1] == 2, "StateC should be active");

        p.process_event(Event<by_name("some_event")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 1, "StateB should be active");
        BOOST_CHECK_MESSAGE(p.current_state()[1] == 2, "StateC should be active");

        p.process_event(Event<by_name("some_event")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 0, "StateA should be active");
        BOOST_CHECK_MESSAGE(p.current_state()[1] == 2, "StateC should be active");

        // force termination
        p.process_event(Event<by_name("terminate_event")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 0, "StateA should be active");
        BOOST_CHECK_MESSAGE(p.current_state()[1] == 3, "TerminalState should be active");

        // no more event processing => no state changes
        p.process_event(Event<by_name("some_event")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 0, "StateA should be active");
        BOOST_CHECK_MESSAGE(p.current_state()[1] == 3, "TerminalState should be active");

    }
}

