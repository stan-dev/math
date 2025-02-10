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
    struct player_ : public msm::front::state_machine_def<player_>
    {
        unsigned int start_playback_counter=0;
        unsigned int can_close_drawer_counter=0;
        unsigned int test_fct_counter=0;
        unsigned int entry_counter = 0;
        unsigned int exit_counter = 0;

        BOOST_MSM_PUML_DECLARE_TABLE(
            R"( 
            @startuml Player
            skinparam linetype polyline
            state Player{
                [*]-> Empty_
                Stopped_     -> Playing_   : play          / TestFct,start_playback    [DummyGuard]
                Stopped_     -> Open_      : open_close_   / open_drawer
                Stopped_     -> Stopped_   : stop_

                Open_        -> Empty_     : open_close_    / close_drawer              [can_close_drawer]
                Empty_       --> Open_     : open_close_    / open_drawer
                Empty_       ---> Stopped_ : cd_detected    / store_cd_info2            [good_disk_format && always_true]
                Playing_     --> Stopped_  : stop_          / stop_playback
                Playing_     -> Paused_    : pause_         / pause_playback
                Playing_     --> Open_     : open_close_    / stop_and_open
                Paused_      -> Playing_   : end_pause      / resume_playback2
                Paused_      --> Stopped_  : stop_          / stop_playback
                Paused_      --> Open_     : open_close_    / stop_and_open
                Playing_ : flag CDLoaded
                Playing_ : entry increment_entry [always_true]
                Playing_ : exit increment_exit [always_true]
                Paused_  : flag CDLoaded
                Stopped_ : flag CDLoaded                
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
    typedef msm::back11::state_machine<player_> player;


    BOOST_AUTO_TEST_CASE(back11_only_string_puml_test)
    {     
        player p;
        static_assert(msm::back11::get_number_of_regions<typename player::initial_state>::type::value == 1);
        static_assert(::boost::mpl::size<typename player::transition_table>::type::value == 12);

        p.start(); 

        p.process_event(Event<by_name("open_close_")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 1,"Open should be active"); //Open
        BOOST_CHECK_MESSAGE(p.is_flag_active<Flag<by_name("CDLoaded")>>() == false, "CDLoaded should not be active");

        p.process_event(Event<by_name("open_close_")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 2,"Empty should be active"); //Empty
        BOOST_CHECK_MESSAGE(p.can_close_drawer_counter == 1,"guard not called correctly");

        p.process_event(Event<by_name("cd_detected")>{"louie, louie", DISK_DVD});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 2,"Empty should be active"); //Empty

        BOOST_CHECK_MESSAGE(p.entry_counter == 0, "Playing entry not called correctly");
        BOOST_CHECK_MESSAGE(p.exit_counter == 0, "Playing exit not called correctly");

        p.process_event(Event<by_name("cd_detected")>{"louie, louie", DISK_CD});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 3,"Playing should be active"); //Playing
        BOOST_CHECK_MESSAGE(p.start_playback_counter == 1,"action not called correctly");
        BOOST_CHECK_MESSAGE(p.test_fct_counter == 1,"action not called correctly");
        BOOST_CHECK_MESSAGE(p.is_flag_active < Flag<by_name("CDLoaded")>>() == true, "CDLoaded should be active");
        BOOST_CHECK_MESSAGE(p.entry_counter == 1, "Playing entry not called correctly");
        BOOST_CHECK_MESSAGE(p.exit_counter == 0, "Playing exit not called correctly");

        p.process_event(Event<by_name("pause_")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 4,"Paused should be active"); //Paused
        BOOST_CHECK_MESSAGE(p.is_flag_active < Flag<by_name("CDLoaded")>>() == true, "CDLoaded should be active");
        BOOST_CHECK_MESSAGE(p.entry_counter == 1, "Playing entry not called correctly");
        BOOST_CHECK_MESSAGE(p.exit_counter == 1, "Playing exit not called correctly");

        // go back to Playing
        p.process_event(Event<by_name("end_pause")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 3,"Playing should be active"); //Playing
        BOOST_CHECK_MESSAGE(p.entry_counter == 2, "Playing entry not called correctly");
        BOOST_CHECK_MESSAGE(p.exit_counter == 1, "Playing exit not called correctly");

        p.process_event(Event<by_name("pause_")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 4,"Paused should be active"); //Paused
        BOOST_CHECK_MESSAGE(p.entry_counter == 2, "Playing entry not called correctly");
        BOOST_CHECK_MESSAGE(p.exit_counter == 2, "Playing exit not called correctly");

        p.process_event(Event<by_name("stop_")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 0,"Stopped should be active"); //Stopped
        BOOST_CHECK_MESSAGE(p.is_flag_active < Flag<by_name("CDLoaded")>>() == true, "CDLoaded should be active");

        p.process_event(Event<by_name("stop_")>{});
        BOOST_CHECK_MESSAGE(p.current_state()[0] == 0,"Stopped should be active"); //Stopped

    }
}

