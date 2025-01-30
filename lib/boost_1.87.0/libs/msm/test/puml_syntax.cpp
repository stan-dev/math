// Copyright 2024 Christophe Henry
// henry UNDERSCORE christophe AT hotmail DOT com
// This is an extended version of the state machine available in the boost::mpl library
// Distributed under the same license as the original.
// Copyright for the original version:
// Copyright 2005 David Abrahams and Aleksey Gurtovoy. Distributed
// under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/ignore_unused.hpp>

#include <boost/msm/front/puml/puml.hpp>

#ifndef BOOST_MSM_NONSTANDALONE_TEST
#define BOOST_TEST_MODULE puml_test
#endif
#include <boost/test/unit_test.hpp>

using namespace boost::msm::front;
using namespace boost::msm::front::puml;

namespace
{

    BOOST_AUTO_TEST_CASE( puml_test )
    {        
        constexpr auto stt_ = R"([*]-> StateA
                            StateA -> StateB: evt1/ A2 [G1]
                            StateB -> StateB : evt1/ A1,A2 [G1]
                            
                            StateB -> StateA: evt1)";
        auto stt =
            boost::msm::front::puml::create_transition_table([&]() {return stt_; });
        static_assert(std::is_same_v<
            decltype(stt),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Guard <by_name("G1")> >,
            Row<State <by_name("StateB")>, Event <by_name("evt1")>, State <by_name("StateB")>, ActionSequence_<boost::fusion::vector< Action<by_name("A1")>, Action<by_name("A2")>>>, Guard <by_name("G1")> >,
            Row<State <by_name("StateB")>, Event <by_name("evt1")>, State <by_name("StateA")>, none, none >
            >
        >);

        constexpr auto stt2_ = R"([*]->StateA
                             StateA -> StateB: evt1/ A2 [G1&&G2]
                             StateB -> StateB : evt1/ A1, defer [G1]
                             
                             StateB -> StateA: evt1)";
        auto stt2 =
            boost::msm::front::puml::create_transition_table([&]() {return stt2_; });
        static_assert(std::is_same_v<
            decltype(stt2),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, And_<Guard <by_name("G1")>, Guard <by_name("G2")>> >,
            Row<State <by_name("StateB")>, Event <by_name("evt1")>, State <by_name("StateB")>, ActionSequence_<boost::fusion::vector< Action<by_name("A1")>, Defer>>, Guard <by_name("G1")> >,
            Row<State <by_name("StateB")>, Event <by_name("evt1")>, State <by_name("StateA")>, none, none >
            >
        >);

        constexpr auto stt3_ = R"(StateA -> StateB: evt1/ A2 [G1 && G2 && G3])";
        auto stt3 =
            boost::msm::front::puml::create_transition_table([&]() {return stt3_; });
        static_assert(std::is_same_v<
            decltype(stt3),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, And_<Guard <by_name("G1")>, And_<Guard <by_name("G2")>, Guard <by_name("G3")>>> >
            >
        >);

        constexpr auto stt4_ = R"(StateA -> StateB: evt1/ A2 [G1 || G2])";
        auto stt4 =
            boost::msm::front::puml::create_transition_table([&]() {return stt4_; });
        static_assert(std::is_same_v<
            decltype(stt4),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Or_<Guard <by_name("G1")>, Guard <by_name("G2")>> >
            >
        >);

        constexpr auto stt5_ = R"(StateA -> StateB: evt1/ A2 [G1 && G2 || G3])";
        auto stt5 =
            boost::msm::front::puml::create_transition_table([&]() {return stt5_; });
        static_assert(std::is_same_v<
            decltype(stt5),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Or_<And_<Guard <by_name("G1")>, Guard <by_name("G2")>>, Guard <by_name("G3")>> >
            >
        >);

        constexpr auto stt6_ = R"(StateA ---> StateB: evt1/ A2 [!G1])";
        auto stt6 =
            boost::msm::front::puml::create_transition_table([&]() {return stt6_; });
        static_assert(std::is_same_v<
            decltype(stt6),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Not_<Guard <by_name("G1")>> >
            >
        >);

        constexpr auto stt7_ = R"(StateA -> StateB: evt1/ A2 [G1 || !G2])";
        auto stt7 =
            boost::msm::front::puml::create_transition_table([&]() {return stt7_; });
        static_assert(std::is_same_v<
            decltype(stt7),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Or_<Guard <by_name("G1")>, Not_<Guard <by_name("G2")>>> >
            >
        >);

        constexpr auto stt8_ = R"(StateA -> StateB: evt1/ A2 [Or(G1,Not(G2))])";
        auto stt8 =
            boost::msm::front::puml::create_transition_table([&]() {return stt8_; });
        static_assert(std::is_same_v<
            decltype(stt8),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Or_<Guard <by_name("G1")>, Not_<Guard <by_name("G2")>>> >
            >
        >);

        auto stt9 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateB: evt1/ A2 [(G1 || G2) && G3])"; });
        static_assert(std::is_same_v<
            decltype(stt9),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, And_<Or_<Guard <by_name("G1")>, Guard <by_name("G2")>>, Guard <by_name("G3")>> >
            >
        >);

        auto stt10 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateB: evt1/ A2 [(G1 && G2) || G3])"; });
        static_assert(std::is_same_v<
            decltype(stt10),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Or_<And_<Guard <by_name("G1")>, Guard <by_name("G2")>>, Guard <by_name("G3")>> >
            >
        >);

        auto stt11 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateB: evt1/ A2 [G3 && ( G1 || G2 )])"; });
        static_assert(std::is_same_v<
            decltype(stt11),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, And_<Guard <by_name("G3")>, Or_<Guard <by_name("G1")>, Guard <by_name("G2")>>> >
            >
        >);

        auto stt12 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateB: evt1/ A2 [G3 || ( G1 && G2 )])"; });
        static_assert(std::is_same_v<
            decltype(stt12),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Or_<Guard <by_name("G3")>, And_<Guard <by_name("G1")>, Guard <by_name("G2")>>> >
            >
        >);

        auto stt13 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateB: evt1/ A2 [!( G1 && G2 )])"; });
        static_assert(std::is_same_v<
            decltype(stt13),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Not_<And_<Guard <by_name("G1")>, Guard <by_name("G2")>>> >
            >
        >);

        auto stt14 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateB: evt1/ A2 [G3 || !( G1 && G2 )])"; });
        static_assert(std::is_same_v<
            decltype(stt14),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>, Or_<Guard <by_name("G3")>, Not_<And_<Guard <by_name("G1")>, Guard <by_name("G2")>>>> >
            >
        >);

        auto stt15 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateB: evt1/ A2 [(G3 && G1) || ( G1 && G2 )])"; });
        static_assert(std::is_same_v<
            decltype(stt15),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>,
            Or_<And_<Guard <by_name("G3")>, Guard <by_name("G1")>>, And_<Guard <by_name("G1")>, Guard <by_name("G2")>>> >
            >
        >);

        auto stt16 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateB: evt1/ A2 [(G3 && G1) || !(G1 && G2)])"; });
        static_assert(std::is_same_v<
            decltype(stt16),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateB")>, Action <by_name("A2")>,
            Or_<And_<Guard <by_name("G3")>, Guard <by_name("G1")>>, Not_<And_<Guard <by_name("G1")>, Guard <by_name("G2")>>>> >
            >
        >);

        // internal detail::Transition
        auto stt17 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateA: -evt1/ A2 [G1])"; });
        static_assert(std::is_same_v<
            decltype(stt17),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, none, Action <by_name("A2")>, Guard <by_name("G1")>>
            >
        >);

        // kleene
        auto stt18 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateA: * / A2 [G1])"; });
        static_assert(std::is_same_v<
            decltype(stt18),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, boost::any, State <by_name("StateA")>, Action <by_name("A2")>, Guard <by_name("G1")>>
            >
        >);

        // internal + kleene
        auto stt19 =
            boost::msm::front::puml::create_transition_table([]() {return R"(StateA -> StateA: -* / A2 [G1])"; });
        static_assert(std::is_same_v<
            decltype(stt19),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, boost::any, none, Action <by_name("A2")>, Guard <by_name("G1")>>
            >
        >);

        // test init states
        static_assert(std::is_same_v<
            decltype(boost::msm::front::puml::create_transition_table([]() {return R"([*] -> StateA
                                                         StateA --> StateA: evt1 / defer)"; })),
            boost::fusion::vector<
            Row<State <by_name("StateA")>, Event <by_name("evt1")>, State <by_name("StateA")>, Defer, none >
            >
        > );


        State< by_name("StateXXX")> s1;
        boost::ignore_unused(s1);

        // initial states
        constexpr auto stt9_ = R"(
                [*] --> StateA
                StateA -> StateB: evt1 [G1]
                StateB -> StateB : evt1 / A1 [G1]
                --
                [*]->StateC
                StateC -> StateD: -*)";
        auto inits =
            create_initial_states([&]() {return stt9_; });
        static_assert(std::is_same_v<
            decltype(inits),
            boost::fusion::vector< State <by_name("StateA")>, State <by_name("StateC")>>
        >);

        //constexpr auto name_ = "StateA";
        //std::cout << typeid(boost::msm::front::puml::detail::parse_flags(stt_, name_)).name() << std::endl;
    }
}

