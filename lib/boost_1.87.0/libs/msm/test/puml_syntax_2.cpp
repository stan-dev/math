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

#ifndef BOOST_MSM_NONSTANDALONE_TEST
#define BOOST_TEST_MODULE puml_test_2
#endif
#include <boost/test/unit_test.hpp>

using namespace boost::msm::front;
using namespace boost::msm::front::puml;

namespace
{

    BOOST_AUTO_TEST_CASE( puml_test_2 )
    {        
        constexpr auto stt20_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]
                            StateA : flag F1
                            StateA -> StateC: evt2/ A2 [G1])";
        auto stt20 =
            boost::msm::front::puml::create_transition_table([&]() {return stt20_; });
        static_assert(std::is_same_v<
            decltype(stt20),
            boost::fusion::vector<
            Row<State <by_name("StateA"), boost::fusion::vector<Flag<by_name("F1")>>>,
                Event <by_name("evt1")>, 
                State <by_name("StateB")>, 
                Action <by_name("A2")>, 
                Guard <by_name("G1")> >,
            Row<State <by_name("StateA"), boost::fusion::vector<Flag<by_name("F1")>>>, 
                Event <by_name("evt2")>, 
                State <by_name("StateC")>, 
                Action <by_name("A2")>, 
                Guard <by_name("G1")> >
            >
        >);

        constexpr auto stt21_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]
                            StateA : flag F1
                            StateA : flag F2)";
        auto stt21 =
            boost::msm::front::puml::create_transition_table([&]() {return stt21_; });
        static_assert(std::is_same_v<
            decltype(stt21),
            boost::fusion::vector<
            Row<State <by_name("StateA"), boost::fusion::vector<Flag<by_name("F1")>, Flag<by_name("F2")>>>, 
                Event <by_name("evt1")>, 
                State <by_name("StateB")>,
                Action <by_name("A2")>, 
                Guard <by_name("G1")> >
            >
        >);

        constexpr auto stt22_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]
                            StateA : flag F1
                            StateB : flag F2)";
        auto stt22 =
            boost::msm::front::puml::create_transition_table([&]() {return stt22_; });
        static_assert(std::is_same_v<
            decltype(stt22),
            boost::fusion::vector<
            Row<State <by_name("StateA"), boost::fusion::vector<Flag<by_name("F1")>>>, 
                Event <by_name("evt1")>, 
                State <by_name("StateB"), boost::fusion::vector<Flag<by_name("F2")>>>,
                Action <by_name("A2")>, 
                Guard <by_name("G1")> >
            >
        >);

       constexpr auto stt23_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]
                            StateA : flag F1
                            StateA : flag F2
                            StateB : flag F2)";
        auto stt23 =
            boost::msm::front::puml::create_transition_table([&]() {return stt23_; });
        static_assert(std::is_same_v<
            decltype(stt23),
            boost::fusion::vector<
            Row<State <by_name("StateA"), boost::fusion::vector<Flag<by_name("F1")>, Flag<by_name("F2")>>>, 
                Event <by_name("evt1")>, 
                State <by_name("StateB"), boost::fusion::vector<Flag<by_name("F2")>>>,
                Action <by_name("A2")>, 
                Guard <by_name("G1")> >
            >
        >);

        constexpr auto stt24_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]
                            StateA : entry A1
                            StateA : entry A2
                            StateB : entry A1
        )";
        auto stt24 =
            boost::msm::front::puml::create_transition_table([&]() {return stt24_; });

        static_assert(std::is_same_v<
            decltype(stt24),
            boost::fusion::vector<
                Row<
                    State < by_name("StateA"), 
                        boost::fusion::vector<>, 
                        boost::fusion::vector<
                            boost::msm::front::puml::detail::pair_type <
                                Action<by_name("A1")>,
                                none
                            >,
                            boost::msm::front::puml::detail::pair_type <
                                Action<by_name("A2")>,
                                none
                            >
                        >
                    >,
                    Event <by_name("evt1")>, 
                    State < by_name("StateB"), 
                            boost::fusion::vector<>, 
                            boost::fusion::vector<
                                boost::msm::front::puml::detail::pair_type <
                                    Action<by_name("A1")>,
                                    none
                                >
                            >>,
                    Action <by_name("A2")>, 
                    Guard <by_name("G1")> >
                >
            >);

        constexpr auto stt25_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]
                            StateA : entry A1 [G1]
                            StateA : entry A2                            
                            StateB : entry A1 [G1]
                            StateB : entry A2
        )";
        auto stt25 =
            boost::msm::front::puml::create_transition_table([&]() {return stt25_; });

        static_assert(std::is_same_v<
            decltype(stt25),
            boost::fusion::vector<
                Row<
                    State < by_name("StateA"), 
                        boost::fusion::vector<>, 
                        boost::fusion::vector<
                            boost::msm::front::puml::detail::pair_type <
                                Action<by_name("A1")>,
                                Guard<by_name("G1")>
                            >,
                            boost::msm::front::puml::detail::pair_type <
                                Action<by_name("A2")>,
                                none
                            >
                        >
                    >,
                    Event <by_name("evt1")>, 
                    State < by_name("StateB"), 
                        boost::fusion::vector<>, 
                        boost::fusion::vector<
                            boost::msm::front::puml::detail::pair_type <
                                Action<by_name("A1")>,
                                Guard<by_name("G1")>
                            >,
                            boost::msm::front::puml::detail::pair_type <
                                Action<by_name("A2")>,
            none
                            >
                        >
                    >,
                    Action <by_name("A2")>, 
                    Guard <by_name("G1")> >
                >
            >);

        constexpr auto stt26_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]
                            StateA : entry A1,A2 [G1 && G2]
                            StateB : entry A1,A2 [G1&&G2])";
        auto stt26 =
            boost::msm::front::puml::create_transition_table([&]() {return stt26_; });

        static_assert(std::is_same_v<
            decltype(stt26),
            boost::fusion::vector<
                Row<
                    State < by_name("StateA"), 
                        boost::fusion::vector<>, 
                        boost::fusion::vector<
                            boost::msm::front::puml::detail::pair_type <
                                boost::msm::front::ActionSequence_< boost::fusion::vector<
                                    Action<by_name("A1")>, Action<by_name("A2")>
                                >>,
                                And_<Guard <by_name("G1")>, Guard <by_name("G2")>>
                            >
                        >
                    >,
                    Event <by_name("evt1")>, 
                    State < by_name("StateB"), 
                        boost::fusion::vector<>, 
                        boost::fusion::vector<
                            boost::msm::front::puml::detail::pair_type <
                                boost::msm::front::ActionSequence_< boost::fusion::vector<
                                    Action<by_name("A1")>, Action<by_name("A2")>
                                >>,
                                And_<Guard <by_name("G1")>, Guard <by_name("G2")>>
                            >
                        >
                    >,
                    Action <by_name("A2")>, 
                    Guard <by_name("G1")> >
                >
            >);

        constexpr auto stt27_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]
                            StateA : exit A1,A2 [G1 && G2]
                            StateB : entry A1,A2 [G1&&G2])";
        auto stt27 =
            boost::msm::front::puml::create_transition_table([&]() {return stt27_; });

        static_assert(std::is_same_v<
            decltype(stt27),
            boost::fusion::vector<
                Row<
                    State < by_name("StateA"), 
                        boost::fusion::vector<>,
                        boost::fusion::vector<>,
                        boost::fusion::vector<
                            boost::msm::front::puml::detail::pair_type <
                                boost::msm::front::ActionSequence_< boost::fusion::vector<
                                    Action<by_name("A1")>, Action<by_name("A2")>
                                >>,
                                And_<Guard <by_name("G1")>, Guard <by_name("G2")>>
                            >
                        >
                    >,
                    Event <by_name("evt1")>, 
                    State < by_name("StateB"), 
                        boost::fusion::vector<>, 
                        boost::fusion::vector<
                            boost::msm::front::puml::detail::pair_type <
                                boost::msm::front::ActionSequence_< boost::fusion::vector<
                                    Action<by_name("A1")>, Action<by_name("A2")>
                                >>,
                                And_<Guard <by_name("G1")>, Guard <by_name("G2")>>
                            >
                        >
                    >,
                    Action <by_name("A2")>, 
                    Guard <by_name("G1")> >
                >
            >);

        constexpr auto stt28_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]                            
                            StateB -> [*])";
        auto stt28 =
            boost::msm::front::puml::create_transition_table([&]() {return stt28_; });
        static_assert(std::is_same_v<
            decltype(stt28),
            boost::fusion::vector<
            Row<State <by_name("StateA")>,
                Event <by_name("evt1")>, 
                State <by_name("StateB"), boost::fusion::vector<boost::msm::TerminateFlag>>,
                Action <by_name("A2")>, 
                Guard <by_name("G1")> >
            >
        >);

        
        constexpr auto stt29_ = R"(   
            @startuml Player
            skinparam linetype polyline
            state Player{         
                [*]-> StateA
                StateA -> StateB: evt1/ A2 [G1]  
                --
                [*]-> StateC
                StateC -> TerminalState : terminate_event                         
                TerminalState -> [*]
            }
            @enduml
        )";
        auto stt29 =
            boost::msm::front::puml::create_transition_table([&]() {return stt29_; });
        static_assert(std::is_same_v<
            decltype(stt29),
            boost::fusion::vector<
            Row<State <by_name("StateA")>,
                Event <by_name("evt1")>, 
                State <by_name("StateB")>,
                Action <by_name("A2")>, 
                Guard <by_name("G1")> >,
            Row<State <by_name("StateC")>,
                Event <by_name("terminate_event")>, 
                State <by_name("TerminalState"), boost::fusion::vector<boost::msm::TerminateFlag>>,
                none, 
                none>
            >
        >);
        static_assert(boost::msm::front::puml::detail::count_terminates(stt29_) == 1);
        static_assert(boost::msm::front::puml::detail::count_inits(stt29_) == 2);

        //std::cout << "dbg1:" << typeid(stt25).name() << std::endl << std::endl;
        //std::cout << "dbg2:" << typeid(stt29).name() << std::endl << std::endl;
    }
}

