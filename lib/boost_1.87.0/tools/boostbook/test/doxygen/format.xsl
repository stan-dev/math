<?xml version="1.0" encoding="utf-8"?>
<!--
   Copyright (c) 2024 Andrey Semashev

   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
-->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="xml" encoding="utf-8"/>
  <xsl:param name="indent-step" select="'  '"/>

  <xsl:template match="processing-instruction() | comment()">
    <xsl:param name="indent" select="''"/>

    <xsl:text>&#10;</xsl:text><!-- Newline -->
    <xsl:value-of select="$indent"/>
    <xsl:copy/>
  </xsl:template>

  <xsl:template match="text()">
    <xsl:param name="indent" select="''"/>

    <xsl:if test="normalize-space(.) != ''">
      <xsl:text>&#10;</xsl:text><!-- Newline -->
      <xsl:value-of select="$indent"/>
      <xsl:value-of select="normalize-space(.)"/>
    </xsl:if>
  </xsl:template>

  <xsl:template match="*">
    <xsl:param name="indent" select="''"/>

    <xsl:text>&#10;</xsl:text><!-- Newline -->
    <xsl:value-of select="$indent"/>
    <xsl:choose>
      <xsl:when test="count(./*) &gt; 0">
        <xsl:copy>
          <xsl:copy-of select="@*"/>
          <xsl:apply-templates>
            <xsl:with-param name="indent" select="concat($indent, $indent-step)"/>
          </xsl:apply-templates>
          <xsl:text>&#10;</xsl:text><!-- Newline -->
          <xsl:value-of select="$indent"/>
        </xsl:copy>
      </xsl:when>
      <xsl:otherwise>
        <xsl:copy-of select="."/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>
</xsl:stylesheet>
