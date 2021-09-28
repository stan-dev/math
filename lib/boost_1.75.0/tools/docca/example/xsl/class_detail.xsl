<!-- CLASS_DETAIL_TEMPLATE BEGIN -->
<xsl:when test="$normal-tparam = 'Allocator'"><xsl:text>__Allocator__</xsl:text></xsl:when>
<xsl:when test="$normal-tparam = 'InputIterator'"><xsl:text>__InputIterator__</xsl:text></xsl:when>
<xsl:when test="$normal-tparam = 'ConstBufferSequence'"><xsl:text>__ConstBufferSequence__</xsl:text></xsl:when>
<xsl:when test="$normal-tparam = 'MutableBufferSequence'"><xsl:text>__MutableBufferSequence__</xsl:text></xsl:when>
<!-- CLASS_DETAIL_TEMPLATE END -->
