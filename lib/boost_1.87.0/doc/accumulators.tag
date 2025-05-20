<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.5" doxygen_gitid="2f6875a5ca481a69a6f32650c77a667f87d25e88">
  <compound kind="struct">
    <name>boost::accumulators::detail::accumulator_set_result</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1accumulator__set__result.html</filename>
    <templarg>typename AccumulatorSet</templarg>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::accumulator_wrapper</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1accumulator__wrapper.html</filename>
    <templarg>typename Accumulator</templarg>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::contains_feature_of_::apply</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1contains__feature__of___1_1apply.html</filename>
    <templarg>typename Accumulator</templarg>
    <base>boost::accumulators::detail::contains_feature_of</base>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::matches_feature::apply</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1matches__feature_1_1apply.html</filename>
    <templarg>typename Accumulator</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::argument_pack_result</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1argument__pack__result.html</filename>
    <templarg>typename Args</templarg>
    <templarg>typename Feature</templarg>
    <base>accumulator_set_result&lt; boost::remove_reference&lt; parameter::binding&lt; boost::remove_const&lt; boost::remove_reference&lt; Args &gt;::type &gt;::type, tag::accumulator &gt;::type &gt;::type, Feature &gt;</base>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::as_feature</name>
    <filename>structboost_1_1accumulators_1_1as__feature.html</filename>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::as_feature_list</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1as__feature__list.html</filename>
    <templarg>typename Features</templarg>
    <templarg>typename Weight</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::as_feature_list&lt; Features, void &gt;</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1as__feature__list_3_01Features_00_01void_01_4.html</filename>
    <templarg>typename Features</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::as_weighted_feature</name>
    <filename>structboost_1_1accumulators_1_1as__weighted__feature.html</filename>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::build_acc_list</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1build__acc__list.html</filename>
    <templarg>typename First</templarg>
    <templarg>typename Last</templarg>
    <templarg>bool is_empty</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::build_acc_list&lt; First, Last, false &gt;</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1build__acc__list_3_01First_00_01Last_00_01false_01_4.html</filename>
    <templarg>typename First</templarg>
    <templarg>typename Last</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::build_acc_list&lt; First, Last, true &gt;</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1build__acc__list_3_01First_00_01Last_00_01true_01_4.html</filename>
    <templarg>typename First</templarg>
    <templarg>typename Last</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::checked_as_weighted_feature</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1checked__as__weighted__feature.html</filename>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::collect_abstract_features</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1collect__abstract__features.html</filename>
    <templarg>typename Features</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::contains_feature_of</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1contains__feature__of.html</filename>
    <templarg>typename Features</templarg>
    <templarg>typename Accumulator</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::contains_feature_of_</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1contains__feature__of__.html</filename>
    <templarg>typename Features</templarg>
    <class kind="struct">boost::accumulators::detail::contains_feature_of_::apply</class>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::dependencies_of</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1dependencies__of.html</filename>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::depends_on</name>
    <filename>structboost_1_1accumulators_1_1depends__on.html</filename>
    <templarg>typename Feature1</templarg>
    <templarg>typename Feature2</templarg>
    <templarg>...</templarg>
    <base>depends_on_base&lt; mpl::transform&lt; mpl::vector&lt; Feature1, Feature2,... &gt;, as_feature&lt; mpl::_1 &gt; &gt;::type &gt;</base>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::depends_on_base</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1depends__on__base.html</filename>
    <templarg>typename Features</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::extractor</name>
    <filename>structboost_1_1accumulators_1_1extractor.html</filename>
    <templarg>typename Feature</templarg>
    <class kind="struct">boost::accumulators::extractor::result</class>
    <class kind="struct">boost::accumulators::extractor::result&lt; this_type(A1)&gt;</class>
    <member kind="function">
      <type>detail::extractor_result&lt; Arg1, Feature &gt;::type</type>
      <name>operator()</name>
      <anchorfile>structboost_1_1accumulators_1_1extractor.html</anchorfile>
      <anchor>a543291f5a932f0d9b55fb627e62e0362</anchor>
      <arglist>(Arg1 const &amp;arg1) const</arglist>
    </member>
    <member kind="function">
      <type>detail::extractor_result&lt; AccumulatorSet, Feature &gt;::type</type>
      <name>operator()</name>
      <anchorfile>structboost_1_1accumulators_1_1extractor.html</anchorfile>
      <anchor>a4756c17711c4ca273522b860dcf03cfe</anchor>
      <arglist>(AccumulatorSet const &amp;acc, A1 const &amp;a1) const</arglist>
    </member>
    <member kind="function">
      <type>detail::extractor_result&lt; AccumulatorSet, Feature &gt;::type</type>
      <name>operator()</name>
      <anchorfile>structboost_1_1accumulators_1_1extractor.html</anchorfile>
      <anchor>a148b1b51c241313d253856e084d3c701</anchor>
      <arglist>(AccumulatorSet const &amp;acc, A1 const &amp;a1, A2 const &amp;a2,...)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::extractor_result</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1extractor__result.html</filename>
    <templarg>typename A</templarg>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::feature_of</name>
    <filename>structboost_1_1accumulators_1_1feature__of.html</filename>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::feature_tag</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1feature__tag.html</filename>
    <templarg>typename Accumulator</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::insert_dependencies</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1insert__dependencies.html</filename>
    <templarg>typename FeatureMap</templarg>
    <templarg>typename Feature</templarg>
    <templarg>typename Weight</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::insert_feature</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1insert__feature.html</filename>
    <templarg>typename FeatureMap</templarg>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::insert_sequence</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1insert__sequence.html</filename>
    <templarg>typename FeatureMap</templarg>
    <templarg>typename Features</templarg>
    <templarg>typename Weight</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::is_dependent_on</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1is__dependent__on.html</filename>
    <templarg>typename A</templarg>
    <templarg>typename B</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::meta::make_acc_list</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1meta_1_1make__acc__list.html</filename>
    <templarg>typename Sequence</templarg>
    <base>build_acc_list&lt; fusion::result_of::begin&lt; Sequence &gt;::type, fusion::result_of::end&lt; Sequence &gt;::type &gt;</base>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::make_accumulator_tuple</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1make__accumulator__tuple.html</filename>
    <templarg>typename Features</templarg>
    <templarg>typename Sample</templarg>
    <templarg>typename Weight</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::matches_feature</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1matches__feature.html</filename>
    <templarg>typename Feature</templarg>
    <class kind="struct">boost::accumulators::detail::matches_feature::apply</class>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::extractor::result</name>
    <filename>structboost_1_1accumulators_1_1extractor_1_1result.html</filename>
    <templarg>typename F</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::extractor::result&lt; this_type(A1)&gt;</name>
    <filename>structboost_1_1accumulators_1_1extractor_1_1result_3_01this__type_07A1_08_4.html</filename>
    <templarg>typename A1</templarg>
    <base>extractor_result&lt; A1, Feature &gt;</base>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::set_insert_range</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1set__insert__range.html</filename>
    <templarg>typename Set</templarg>
    <templarg>typename Range</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::to_accumulator</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1to__accumulator.html</filename>
    <templarg>typename Feature</templarg>
    <templarg>typename Sample</templarg>
    <templarg>typename Weight</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::to_accumulator&lt; Feature, Sample, tag::external&lt; Weight, Tag, AccumulatorSet &gt; &gt;</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1to__accumulator_3_01Feature_00_01Sample_00_01tag_1_1ext0c26d8c3c18bcca084cb467b003ed836.html</filename>
    <templarg>typename Feature</templarg>
    <templarg>typename Sample</templarg>
    <templarg>typename Weight</templarg>
    <templarg>typename Tag</templarg>
    <templarg>typename AccumulatorSet</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::undroppable</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1undroppable.html</filename>
    <templarg>typename Feature</templarg>
  </compound>
  <compound kind="struct">
    <name>boost::accumulators::detail::undroppable&lt; tag::droppable&lt; Feature &gt; &gt;</name>
    <filename>structboost_1_1accumulators_1_1detail_1_1undroppable_3_01tag_1_1droppable_3_01Feature_01_4_01_4.html</filename>
    <templarg>typename Feature</templarg>
  </compound>
</tagfile>
