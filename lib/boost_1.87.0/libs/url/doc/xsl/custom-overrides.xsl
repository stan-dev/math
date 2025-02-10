<xsl:stylesheet version="3.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema"
  exclude-result-prefixes="xs"
  expand-text="yes">

  <xsl:variable name="doc-ref" select="'url.ref'"/>
  <xsl:variable name="doc-ns" select="'boost::urls'"/>
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

  <xsl:template mode="convenience-header" match="@file[contains(., 'boost/url/grammar')]" priority="1">url/grammar.hpp</xsl:template>
  <xsl:template mode="convenience-header" match="@file[contains(., 'boost/url')]">url.hpp</xsl:template>
  <xsl:template mode="convenience-header" match="@file"/>

  <xsl:variable name="emphasized-template-parameter-types" select="
    'Allocator',
    'CharSet',
    'Rule', 'Rule1', 'Rule2', 'Rules',
    'StringToken'
  "/>

</xsl:stylesheet>
