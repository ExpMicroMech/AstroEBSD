<?xml version="1.0" encoding="utf-8"?>

<!--
This is an XSL stylesheet which converts mscript XML files into DocBook.
Use the XSLT command to perform the conversion.

Copyright 1984-2020 The MathWorks, Inc.
-->

<!DOCTYPE xsl:stylesheet [ <!ENTITY nbsp "&#160;"> <!ENTITY reg "&#174;"> ]>
<xsl:stylesheet
  version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:mwsh="http://www.mathworks.com/namespace/mcode/v1/syntaxhighlight.dtd"
  exclude-result-prefixes="mwsh">
  <xsl:output method="xml"
    indent="no"
    doctype-system="http://www.oasis-open.org/docbook/xml/4.2/docbookx.dtd"
    doctype-public="-//OASIS//DTD DocBook XML V4.2//EN"/>
  <xsl:strip-space elements="mwsh:code"/>

<xsl:variable name="title">
  <xsl:variable name="dTitle" select="//steptitle[@style='document']"/>
  <xsl:choose>
    <xsl:when test="$dTitle"><xsl:value-of select="$dTitle"/></xsl:when>
    <xsl:otherwise><xsl:value-of select="mscript/m-file"/></xsl:otherwise>
  </xsl:choose>
</xsl:variable>


<xsl:template match="mscript">
<article>
  <xsl:attribute name="xml:base">
    <xsl:text>file:///</xsl:text>
    <xsl:call-template name="globalReplace">
      <xsl:with-param name="outputString">
        <xsl:value-of select="outputdir"/>
      </xsl:with-param>
      <xsl:with-param name="target" select="'\'"/>
      <xsl:with-param name="replacement" select="'/'"/>
    </xsl:call-template>
    <xsl:text>/</xsl:text>
  </xsl:attribute>

  <!-- head -->
<xsl:comment>
This HTML is auto-generated from an M-file.
To make changes, update the M-file and republish this document.
</xsl:comment>

    <xsl:call-template name="header"/>

    <!-- Determine if the there should be an introduction section. -->
    <xsl:variable name="hasIntro" select="count(cell[@style = 'overview'])"/>

    <!-- If there is an introduction, display it. -->
    <xsl:if test = "$hasIntro">
            <title><xsl:apply-templates select="cell[1]/steptitle"/></title>
        <abstract><xsl:apply-templates select="cell[1]/text"/></abstract>
    </xsl:if>

    <xsl:variable name="sections" select="cell[not(@style = 'overview') and steptitle] | cell[not(@style = 'overview')][1]"/>

    <!-- Loop over each cell -->

    <xsl:for-each select="$sections">
        <section>
            <title>
                <xsl:apply-templates select="steptitle"/>
            </title>

            <xsl:apply-templates select="text"/>
            <xsl:apply-templates select="mcode-xmlized"/>
            <xsl:apply-templates select="mcodeoutput|img"/>

            <xsl:choose>
                <xsl:when test="following-sibling::cell[steptitle]">
                    <xsl:for-each select="
                                  following-sibling::*[
                                  count(. | current()/following-sibling::cell[steptitle][1]/preceding-sibling::cell) =
                                  count(current()/following-sibling::cell[steptitle][1]/preceding-sibling::cell)
                    ]">
                        <xsl:apply-templates select="text"/>
                        <xsl:apply-templates select="mcode-xmlized"/>
                        <xsl:apply-templates select="mcodeoutput|img"/>
                    </xsl:for-each>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:for-each select="following-sibling::cell">
                        <xsl:apply-templates select="text"/>
                        <xsl:apply-templates select="mcode-xmlized"/>
                        <xsl:apply-templates select="mcodeoutput|img"/>
                    </xsl:for-each>
                </xsl:otherwise>
            </xsl:choose>

        </section>
    </xsl:for-each>

    <xsl:call-template name="footer"/>

</article>
</xsl:template>

<!-- Header -->
<xsl:template name="header">
</xsl:template>

<!-- Footer -->
<xsl:template name="footer">
  <literallayout>
<emphasis><xsl:value-of select="copyright"/>
Published with MATLAB&reg; R<xsl:value-of select="release"/></emphasis>
  </literallayout>
</xsl:template>

<!-- HTML Tags in text sections -->
<xsl:template match="p">
  <para><xsl:apply-templates/></para>
</xsl:template>
<xsl:template match="ul">
  <itemizedlist><xsl:apply-templates/></itemizedlist>
</xsl:template>
<xsl:template match="ol">
  <orderedlist><xsl:apply-templates/></orderedlist>
</xsl:template>
<xsl:template match="li">
  <listitem><xsl:apply-templates/></listitem>
</xsl:template>
<xsl:template match="pre">
    <code><literallayout><xsl:apply-templates/></literallayout></code>
</xsl:template>
<xsl:template match="b">
  <emphasis role="bold"><xsl:apply-templates/></emphasis >
</xsl:template>
<xsl:template match="i">
  <emphasis><xsl:apply-templates/></emphasis>
</xsl:template>
<xsl:template match="tt">
  <literal><xsl:apply-templates/></literal>
</xsl:template>
<xsl:template match="a">
  <ulink>
    <xsl:attribute name="url"><xsl:value-of select="@href"/></xsl:attribute>
    <xsl:apply-templates/>
  </ulink>
</xsl:template>
<xsl:template match="html"/>
<xsl:template match="latex"/>

<!-- Code input and output -->

<xsl:template match="mcode-xmlized">
  <programlisting><xsl:apply-templates/><xsl:text><!-- g162495 -->
</xsl:text></programlisting>
</xsl:template>

<xsl:template match="mcodeoutput">
  <xsl:choose>
    <xsl:when test="substring(.,0,7)='&lt;html&gt;'">
    </xsl:when>
    <xsl:when test="./@class = 'codeoutput error'">
        <!-- This is an error; so we color it red. -->
        <screen><computeroutput><emphasis role="red"><xsl:apply-templates/></emphasis></computeroutput></screen>
    </xsl:when>
    <xsl:otherwise>
        <!-- Ordinary MATLAB output, not an error -->
        <screen><computeroutput><emphasis><xsl:apply-templates/></emphasis></computeroutput></screen>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<!-- Figure and model snapshots and equations -->
<xsl:template match="img[@class='equation']">
  <inlineequation>
    <inlinemediaobject>
      <imageobject>
        <imagedata scalefit="1">
          <xsl:variable name="space_escaped">
            <xsl:call-template name="globalReplace">
              <xsl:with-param name="outputString">
                <xsl:value-of select="@src"/>
              </xsl:with-param>
              <xsl:with-param name="target" select="' '"/>
              <xsl:with-param name="replacement" select="'%20'"/>
            </xsl:call-template>
          </xsl:variable>

          <xsl:attribute name="fileref">
            <xsl:call-template name="globalReplace">
              <xsl:with-param name="outputString">
                <xsl:value-of select="concat('file:///', $space_escaped)"/>
              </xsl:with-param>
              <xsl:with-param name="target" select="'\'"/>
              <xsl:with-param name="replacement" select="'/'"/>
            </xsl:call-template>
          </xsl:attribute>
          <xsl:attribute name="depth"><xsl:value-of select="@height"/></xsl:attribute>
        </imagedata>
      </imageobject>
    </inlinemediaobject>
  </inlineequation>
</xsl:template>

<xsl:template match="img">
  <mediaobject>
    <imageobject>
      <xsl:variable name="width">
        <xsl:choose>
          <xsl:when test="@width">
            <xsl:variable name="in" select="substring-before(@width,'px') div 72"/>
            <xsl:choose>
                <xsl:when test="$in > 5.8333">
                <xsl:text>5.8333</xsl:text>
              </xsl:when>
              <xsl:otherwise>
                  <xsl:value-of select="$in"/>
              </xsl:otherwise>
            </xsl:choose>
            <xsl:text>in</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>100%</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <imagedata scalefit="1" width="{$width}" contentdepth="100%">
        <xsl:variable name="space_escaped">
          <xsl:call-template name="globalReplace">
            <xsl:with-param name="outputString">
              <xsl:value-of select="@src"/>
            </xsl:with-param>
            <xsl:with-param name="target" select="' '"/>
            <xsl:with-param name="replacement" select="'%20'"/>
          </xsl:call-template>
        </xsl:variable>

        <xsl:attribute name="fileref">
          <xsl:call-template name="globalReplace">
            <xsl:with-param name="outputString">
              <xsl:value-of select="concat('file:///', $space_escaped)"/>
            </xsl:with-param>
            <xsl:with-param name="target" select="'\'"/>
            <xsl:with-param name="replacement" select="'/'"/>
           </xsl:call-template>
        </xsl:attribute>
      </imagedata>
    </imageobject>
  </mediaobject>
</xsl:template>

<!-- Stash original code in HTML for easy slurping later. -->

<xsl:template match="originalCode">
  <xsl:variable name="xcomment">
    <xsl:call-template name="globalReplace">
      <xsl:with-param name="outputString" select="."/>
      <xsl:with-param name="target" select="'--'"/>
      <xsl:with-param name="replacement" select="'REPLACE_WITH_DASH_DASH'"/>
    </xsl:call-template>
  </xsl:variable>
<xsl:comment>
##### SOURCE BEGIN #####
<xsl:value-of select="$xcomment"/>
##### SOURCE END #####
</xsl:comment>
</xsl:template>

<!-- Colors for syntax-highlighted input code -->

<xsl:template match="mwsh:code">
  <xsl:apply-templates/>
</xsl:template>
<xsl:template match="mwsh:keywords">
  <mwsh:keywords><xsl:value-of select="."/></mwsh:keywords>
</xsl:template>
<xsl:template match="mwsh:strings">
  <mwsh:strings><xsl:value-of select="."/></mwsh:strings>
</xsl:template>
<xsl:template match="mwsh:comments">
  <mwsh:comments><xsl:value-of select="."/></mwsh:comments>
</xsl:template>
<xsl:template match="mwsh:unterminated_strings">
  <mwsh:unterminated_strings><xsl:value-of select="."/></mwsh:unterminated_strings>
</xsl:template>
<xsl:template match="mwsh:system_commands">
  <mwsh:system_commands><xsl:value-of select="."/></mwsh:system_commands>
</xsl:template>
<xsl:template match="mwsh:type_section">
  <mwsh:type_section><xsl:value-of select="."/></mwsh:type_section>
</xsl:template>


<!-- Footer information -->

<xsl:template match="copyright">
  <xsl:value-of select="."/>
</xsl:template>
<xsl:template match="revision">
  <xsl:value-of select="."/>
</xsl:template>

<!-- Search and replace  -->
<!-- From http://www.xml.com/lpt/a/2002/06/05/transforming.html -->

<xsl:template name="globalReplace">
  <xsl:param name="outputString"/>
  <xsl:param name="target"/>
  <xsl:param name="replacement"/>
  <xsl:choose>
    <xsl:when test="contains($outputString,$target)">
      <xsl:value-of select=
        "concat(substring-before($outputString,$target),$replacement)"/>
      <xsl:call-template name="globalReplace">
        <xsl:with-param name="outputString"
          select="substring-after($outputString,$target)"/>
        <xsl:with-param name="target" select="$target"/>
        <xsl:with-param name="replacement"
          select="$replacement"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="$outputString"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

</xsl:stylesheet>
