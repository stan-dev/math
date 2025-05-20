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

// back-end
#include <boost/msm/back11/state_machine.hpp>
//front-end
#include <boost/msm/front/state_machine_def.hpp>
#include <boost/msm/front/puml/puml.hpp>

using namespace std;
namespace msm = boost::msm;
using namespace msm::front;
using namespace msm::front::puml;

namespace {
    // Flags. Allow information about a property of the current state
    struct PlayingPaused {};
    struct CDLoaded {};
    struct FirstSongPlaying {};
}

namespace boost::msm::front::puml {
    // events, only if non-default structs
    template<>
    struct Event<by_name("cd_detected")>
    {
        Event(std::string name)
            : name(name)
        {}
        std::string name;
    };
    // The list of FSM states, only if non-default
    template<>
    struct State<by_name("Song1")> : public msm::front::state<>
    {
        typedef boost::fusion::vector<FirstSongPlaying>      flag_list;
        // optional
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { }
    };


    //actions. Must be provided
    template<>
    struct Action<by_name("start_next_song")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("start_prev_song")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    // guard conditions. Must be provided
    template<>
    struct Guard<by_name("start_prev_song_guard")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
            return true;
        }
    };
}

namespace
{
    // front-end: define the FSM structure 
    struct playing_ : public msm::front::state_machine_def<playing_>
    {
        typedef boost::fusion::vector<PlayingPaused, CDLoaded>        flag_list;
        // optional
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { }

        BOOST_MSM_PUML_DECLARE_TABLE(
            R"( 
            @startuml Player
            skinparam linetype polyline
            state Player{
                [*]-> Song1
                Song1     -> Song2   : NextSong         
                Song2     -> Song1   : PreviousSong   / start_prev_song [start_prev_song_guard]
                Song2     -> Song3   : NextSong       / start_next_song
                Song3     -> Song2   : PreviousSong                     [start_prev_song_guard]
            }
            @enduml
        )"
        )

        // Replaces the default no-transition response.
        template <class FSM, class Event>
        void no_transition(Event const&, FSM&, int)
        {
        }
    };
}
namespace boost::msm::front::puml {
    // define submachine with desired back-end
    template<>
    struct State<by_name("PlayingFsm")> : public msm::back11::state_machine<playing_>
    {
    };
}

namespace boost::msm::front::puml {
    // states, only if non-default
    template<>
    struct State<by_name("Paused")> : public msm::front::state<>
    {
        typedef boost::fusion::vector<PlayingPaused, CDLoaded>        flag_list;

        // optional
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { }
    };
    template<>
    struct State<by_name("Stopped")> : public msm::front::state<>
    {
        typedef boost::fusion::vector<CDLoaded>      flag_list;

        // optional
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { }
    };
    template<>
    struct State<by_name("ErrorMode")> : 
        public msm::front::interrupt_state < Event<by_name("end_error")>>
    {
        // optional
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { }
    };
    template<>
    struct State<by_name("ErrorTerminate")> : public msm::front::terminate_state<>
    {
        // optional
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { }
    };

    //actions. Must be provided
    template<>
    struct Action<by_name("open_drawer")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("close_drawer")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("report_error")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("report_end_error")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("my_defer")>
    {
        // mark as deferring to avoid stack overflows in certain conditions
        typedef int deferring_action;
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM& fsm, SourceState&, TargetState&) const
        {
            Event<by_name("play")> e{};
            fsm.defer_event(Event<by_name("play")>{});
        }
    };
    template<>
    struct Action<by_name("start_playback")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("stop_playback")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("pause_playback")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("resume_playback")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT&, FSM&, SourceState&, TargetState&)
        {
        }
    };

    // guard conditions. Must be provided
    template<>
    struct Guard<by_name("is_play_event")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const& evt, FSM& , SourceState&, TargetState&)
        {
            bool is_play = boost::any_cast<Event<by_name("play")>>(&evt) != 0;
            return is_play;
        }
    };
}

namespace
{
    // front-end: define the FSM structure 
    struct player_ : public msm::front::state_machine_def<player_>
    {
        // we want deferred events and no state requires deferred events (only the fsm in the
        // transition table), so the fsm does.
        typedef int activate_deferred_events;

        BOOST_MSM_PUML_DECLARE_TABLE(
            R"( 
            @startuml Player
            skinparam linetype polyline
            state Player{
                [*]-> Empty
                Stopped     -> PlayingFsm     : play          / start_playback
                Stopped     -> Open           : open_close    / open_drawer
                Stopped     -> Stopped        : stop

                Open        -> Empty          : open_close  
                'default defer is provided                          
                Open        -> Open           : -play         / defer
                
                Empty       --> Open          : open_close    / open_drawer
                Empty       ---> Stopped      : cd_detected
                '* is any event (Kleene)
                Empty       -> Empty          : -*            / my_defer                   [is_play_event]                          
                
                PlayingFsm  --> Stopped       : stop          / stop_playback
                PlayingFsm  -> Paused         : pause         / pause_playback
                PlayingFsm  --> Open          : open_close    / stop_playback, open_drawer

                Paused      -> PlayingFsm     : end_pause     / resume_playback
                Paused      --> Stopped       : stop          / stop_playback
                Paused      --> Open          : open_close    / stop_playback, open_drawer
                --
                ' a second region
                [*]-> AllOk
                AllOk       -> ErrorMode      : error_found   / report_error
                ErrorMode   -> AllOk          : end_error     / report_end_error
                AllOk       -> ErrorTerminate : do_terminate
            }
            @enduml
        )"
        )

        // Replaces the default no-transition response.
        template <class FSM,class Event>
        void no_transition(Event const&, FSM&,int)
        {
        }
    };
    // Pick a back-end
    typedef msm::back11::state_machine<player_> player;


    void test( )
    {
        player p;
        // needed to start the highest-level SM. This will call on_entry and mark the start of the SM
        p.start(); 
        // test deferred event
        // deferred in Empty and Open, will be handled only after event cd_detected
        p.process_event(Event<by_name("play")>{});
        assert(p.current_state()[0] == 2); //Empty
        //flags
        assert(p.is_flag_active<CDLoaded>() == false);
        p.process_event(Event<by_name("open_close")>{});
        assert(p.current_state()[0] == 1); //Open
        p.process_event(Event<by_name("open_close")>{});
        assert(p.current_state()[0] == 2); //Empty
        
        //deferred event should have been processed
        p.process_event(Event<by_name("cd_detected")>{"louie, louie"});
        assert(p.current_state()[0] == 3); //Playing
        assert(p.get_state<State<by_name("PlayingFsm")>&>().current_state()[0] == 0);//Song1

        //flags
        assert(p.is_flag_active<PlayingPaused>() == true);
        assert(p.is_flag_active<FirstSongPlaying>() == true);

        p.process_event(Event<by_name("NextSong")>{});
        assert(p.current_state()[0] == 3); //Playing
        assert(p.get_state<State<by_name("PlayingFsm")>&>().current_state()[0] == 1);//Song2

        p.process_event(Event<by_name("NextSong")>{});
        assert(p.current_state()[0] == 3); //Playing
        assert(p.get_state<State<by_name("PlayingFsm")>&>().current_state()[0] == 2);// Song3

        p.process_event(Event<by_name("PreviousSong")>{});
        assert(p.current_state()[0] == 3); //Playing
        assert(p.get_state<State<by_name("PlayingFsm")>&>().current_state()[0] == 1);// Song2
        //flags
        assert(p.is_flag_active<PlayingPaused>() == true);
        assert(p.is_flag_active<FirstSongPlaying>() == false);

        p.process_event(Event<by_name("pause")>{});
        assert(p.current_state()[0] == 4); //Paused
        //flags
        assert(p.is_flag_active<PlayingPaused>() == true);

        // go back to Playing
        p.process_event(Event<by_name("end_pause")>{});
        assert(p.current_state()[0] == 3); //Playing

        p.process_event(Event<by_name("pause")>{});
        assert(p.current_state()[0] == 4); //Paused

        p.process_event(Event<by_name("stop")>{});
        assert(p.current_state()[0] == 0); //Stopped
        //flags
        assert(p.is_flag_active<PlayingPaused>() == false);
        assert(p.is_flag_active<CDLoaded>() == true);

        p.process_event(Event<by_name("stop")>{});
        assert(p.current_state()[0] == 0); //Stopped

        //test interrupt
        assert(p.current_state()[1] == 5); //AllOk
        p.process_event(Event<by_name("error_found")>{});
        assert(p.current_state()[1] == 6); //ErrorMode
        // try generating more events
        p.process_event(Event<by_name("play")>{});
        assert(p.current_state()[1] == 6); //ErrorMode
        assert(p.current_state()[0] == 0); //Stopped
        p.process_event(Event<by_name("end_error")>{});
        assert(p.current_state()[1] == 5); //AllOk
        p.process_event(Event<by_name("play")>{});
        assert(p.current_state()[0] == 3); //Playing

        //test terminate
        assert(p.current_state()[1] == 5); //AllOk
        p.process_event(Event<by_name("do_terminate")>{});
        assert(p.current_state()[1] == 7); //ErrorTerminate
        // try generating more events
        p.process_event(Event<by_name("stop")>{});
        assert(p.current_state()[1] == 7); //ErrorTerminate
        assert(p.current_state()[0] == 3); //Playing

    }
}
int main()
{
    test();
    return 0;
}

