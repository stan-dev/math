// Copyright 2024 Christophe Henry
// henry UNDERSCORE christophe AT hotmail DOT com
// This is an extended version of the state machine available in the boost::mpl library
// Distributed under the same license as the original.
// Copyright for the original version:
// Copyright 2005 David Abrahams and Aleksey Gurtovoy. Distributed
// under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/msm/front/puml/puml.hpp>

using namespace boost::msm::front;
using namespace boost::msm::front::puml;
namespace boost::msm::front::puml {
    // events
    template<>
    struct Event<by_name("play")> { int cpt_ = 0; };
    template<>
    struct Event<by_name("end_pause")> { int cpt_ = 0; };
    template<>
    struct Event<by_name("stop")> {};
    template<>
    struct Event<by_name("pause")> {};
    template<>
    struct Event<by_name("open_close")> {};
    template<>
    struct Event<by_name("internal_evt")> {};
    template<>
    struct Event<by_name("to_ignore")> {};
    // A "complicated" event type that carries some data.
    enum DiskTypeEnum
    {
        DISK_CD = 0,
        DISK_DVD = 1
    };
    template<>
    struct Event<by_name("cd_detected")>
    {
        Event(std::string name, DiskTypeEnum diskType= DISK_CD)
            : name(name),
            disc_type(diskType)
        {}

        std::string name;
        DiskTypeEnum disc_type;
    };

    // The list of FSM states
    template<>
    struct State<by_name("Empty")> : public msm::front::state<>
    {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
        unsigned int empty_internal_guard_counter = 0;
        unsigned int empty_internal_action_counter = 0;
    };
    template<>
    struct State<by_name("Open")> : public msm::front::state<>
    {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
    };
    template<>
    struct State<by_name("Stopped")> : public msm::front::state<>
    {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
    };
    template<>
    struct State<by_name("Paused")> : public msm::front::state<>
    {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
    };
    template<>
    struct State<by_name("Playing")> : public msm::front::state<>
    {
        template <class Event, class FSM>
        void on_entry(Event const&, FSM&) { ++entry_counter; }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
    };
    template<>
    struct State<by_name("Playing2")> : public msm::front::state<>
    {
        template <class Event, class FSM>
        void on_entry(Event const& e, FSM&)
        {
            ++entry_counter;
            event_counter = e.cpt_;
        }
        template <class Event, class FSM>
        void on_exit(Event const&, FSM&) { ++exit_counter; }
        int entry_counter = 0;
        int exit_counter = 0;
        int event_counter = 0;
    };
    //actions
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
    struct Action<by_name("start_playback")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            ++fsm.start_playback_counter;
        }
    };
    template<>
    struct Action<by_name("store_cd_info")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("store_cd_info2")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            fsm.process_event(Event<by_name("play")>{});
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
    template<>
    struct Action<by_name("resume_playback2")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT& e, FSM&, SourceState&, TargetState&)
        {
            ++e.cpt_;
        }
    };
    template<>
    struct Action<by_name("stop_and_open")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("stopped_again")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState&, TargetState&)
        {
        }
    };
    template<>
    struct Action<by_name("internal_action")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            ++fsm.internal_action_counter;
        }
    };
    template<>
    struct Action<by_name("internal_action_fct")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT const&, FSM&, SourceState& src, TargetState&)
        {
            ++src.empty_internal_action_counter;
        }
    };

    template<>
    struct Action<by_name("TestFct")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT& e, FSM& fsm, SourceState&, TargetState&)
        {
            ++e.cpt_;
            ++fsm.test_fct_counter;
        }
    };
    template<>
    struct Action<by_name("increment_entry")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT&, FSM& fsm, SourceState&, TargetState&)
        {
            ++fsm.entry_counter;
        }
    };
    template<>
    struct Action<by_name("increment_exit")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        void operator()(EVT&, FSM& fsm, SourceState&, TargetState&)
        {
            ++fsm.exit_counter;
        }
    };

    // guard conditions
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
    template<>
    struct Guard<by_name("can_close_drawer")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            ++fsm.can_close_drawer_counter;
            return true;
        }
    };
    template<>
    struct Guard<by_name("internal_guard")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            ++fsm.internal_guard_counter;
            return false;
        }
    };
    template<>
    struct Guard<by_name("internal_guard_fct")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM&, SourceState& src, TargetState&)
        {
            ++src.empty_internal_guard_counter;
            return false;
        }
    };
    template<>
    struct Guard<by_name("internal_guard2")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM& fsm, SourceState&, TargetState&)
        {
            ++fsm.internal_guard_counter;
            return true;
        }
    };
    template<>
    struct Guard<by_name("DummyGuard")>
    {
        template <class EVT, class FSM, class SourceState, class TargetState>
        bool operator()(EVT const&, FSM&, SourceState&, TargetState&)
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
}

