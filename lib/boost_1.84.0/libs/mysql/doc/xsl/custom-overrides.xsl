<!--
    Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
   
    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
-->

<xsl:stylesheet version="3.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema" exclude-result-prefixes="xs" expand-text="yes">

  <xsl:variable name="doc-ns" select="'boost::mysql'"/>
  <xsl:variable name="doc-ref" select="'mysql.ref'"/>
  <xsl:variable name="debug" select="0"/>
  <xsl:variable name="private" select="0"/>
  <xsl:variable name="include-private-members" select="false()"/>

  <xsl:template mode="includes-template-footer" match="location">
    <xsl:variable name="convenience-header" as="xs:string?">
      <xsl:apply-templates mode="convenience-header" select="@file"/>
    </xsl:variable>
    <xsl:if test="$convenience-header">
      <xsl:text>{$nl}</xsl:text>
      <xsl:text>Convenience header [include_file boost/{$convenience-header}]</xsl:text>
      <xsl:text>{$nl}</xsl:text>
    </xsl:if>
  </xsl:template>

  <xsl:template mode="convenience-header" match="@file[contains(., 'boost/mysql')]">mysql.hpp</xsl:template>

  <xsl:variable name="emphasized-template-parameter-types" select="
    'CompletionToken',
    'Stream',
    'SocketStream',
    'Executor',
    'WritableFieldTuple',
    'FieldViewFwdIterator',
    'ExecutionRequest',
    'StaticRow',
    'ResultsType',
    'ExecutionStateType'
  "/>

</xsl:stylesheet>
