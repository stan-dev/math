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
// back-end
#include <boost/msm/back11/state_machine.hpp>
//front-end
#include <boost/msm/front/state_machine_def.hpp>
#include <boost/msm/front/puml/puml.hpp>

using namespace std;
namespace msm = boost::msm;
using namespace msm::front;
using namespace msm::front::puml;

// State, events, if non-default must be specialized in boost::msm::front::puml namespace
// Actions and guard must be specialized in boost::msm::front::puml namespace
namespace boost::msm::front::puml {

    // A "complicated" event type that carries some data.
    enum DiskTypeEnum
    {
        DISK_CD = 0,
        DISK_DVD = 1
    };
    template<>
    struct Event<by_name("cd_detected")>
    {
        Event(std::string name, DiskTypeEnum diskType = DISK_CD)
            : name(name),
            disc_type(diskType)
        {}

        std::string name;
        DiskTypeEnum disc_type;
    };

    //Actions
    template<>
    struct Action<by_name("open_drawer")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
            std::cout << "open drawer" << std::endl;
        }
    };
    template<>
    struct Action<by_name("close_drawer")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
            std::cout << "close drawer" << std::endl;
        }
    };
    template<>
    struct Action<by_name("start_playback")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            std::cout << "start playback" << std::endl;
        }
    };
    template<>
    struct Action<by_name("stop_playback")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            std::cout << "stop playback" << std::endl;
        }
    };
    template<>
    struct Action<by_name("resume_playback")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            std::cout << "resume playback" << std::endl;
        }
    };
    template<>
    struct Action<by_name("pause_playback")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            std::cout << "pause playback" << std::endl;
        }
    };
    template<>
    struct Action<by_name("store_cd_info")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            fsm.process_event(Event<by_name("play")>{});
        }
    };
    // Guards
    template<>
    struct Guard<by_name("can_close_drawer")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            return true;
        }
    };
    template<>
    struct Guard<by_name("always_true")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
            return true;
        }
    };
    template<>
    struct Guard<by_name("good_disk_format")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT& evt, FSM&, SourceState&, TargetState&)
        {
            // to test a guard condition, let's say we understand only CDs, not DVD
            if (evt.disc_type != DISK_CD)
            {
                return false;
            }
            return true;
        }
    };
}

namespace
{
    // front-end: define the FSM structure 
    struct player_ : public msm::front::state_machine_def<player_>
    {
        // here is the whole FSM structure defined:
        // Initial states [*]
        // Transitions with actions starting with / and separated by ,
        // and guards between []. Supported are !, &&, || and ()
        // State Entry / Exit with guards
        // Flags
        // -> can have different lengths for nicer PlantUML Viewer preview

        BOOST_MSM_PUML_DECLARE_TABLE(
            R"( 
            @startuml Player
            skinparam linetype polyline
            state Player{
                [*]-> Empty
                Stopped     -> Playing   : play 
                Stopped     -> Open      : open_close    / open_drawer
                Stopped     -> Stopped   : stop

                Open        -> Empty     : open_close    / close_drawer               [can_close_drawer]
                Empty       --> Open     : open_close    / open_drawer
                Empty       ---> Stopped : cd_detected   / store_cd_info              [good_disk_format && always_true]
                Playing     --> Stopped  : stop          / stop_playback
                Playing     -> Paused    : pause
                Playing     --> Open     : open_close    / stop_playback, open_drawer
                Paused      -> Playing   : end_pause     / resume_playback
                Paused      --> Stopped  : stop          / stop_playback
                Paused      --> Open     : open_close    / stop_playback, open_drawer

                Playing : flag CDLoaded
                Playing : entry start_playback [always_true]
                Paused  : entry pause_playback
                Paused  : flag CDLoaded
                Stopped : flag CDLoaded                
            }
            @enduml
        )"
        )

        // Replaces the default no-transition response.
        template <class FSM, class Event>
        void no_transition(Event const& e, FSM&, int state)
        {
            std::cout << "no transition from state " << state
                      << " on event " << typeid(e).name() << std::endl;
        }
    };
    // Pick a back-end
    typedef msm::back11::state_machine<player_> player;


    void test()
    {     
        player p;
        p.start(); 

        p.process_event(Event<by_name("open_close")>{});
        // Open is active

        p.process_event(Event<by_name("open_close")>{});
        // Empty is active

        p.process_event(Event<by_name("cd_detected")>{"louie, louie", DISK_DVD});
        // Empty is active because disk format rejected

   
        p.process_event(Event<by_name("cd_detected")>{"louie, louie", DISK_CD});
        // Playing is active because disk format accepted, flag is also set
        assert(p.is_flag_active < Flag<by_name("CDLoaded")>>());

        
        p.process_event(Event<by_name("pause")>{});
        // Paused is active


        // go back to Playing
        p.process_event(Event<by_name("end_pause")>{});
        // Playing is active

        p.process_event(Event<by_name("stop")>{});
        // Stopped is active

        p.stop();
    }
}
int main()
{
    test();
    return 0;
}


