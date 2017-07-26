##
## Process this file with automake to produce Makefile.in
##

AUTOMAKE_OPTIONS = no-dependencies

DATABASE_FILES = \
  db_general.csv \
  db_other_genomes.csv

SOURCE_CONFIG_FILES = \
  AME_config.xml.in \
  CENTRIMO_config.xml.in \
  DREME_config.xml.in \
  FIMO_config.xml.in \
  GLAM2_config.xml.in \
  GLAM2SCAN_config.xml.in \
  GOMO_config.xml.in \
  MAST_config.xml.in \
  MCAST_config.xml.in \
  MEMECHIP_config.xml.in \
  MEME_config.xml.in \
  MOMO_config.xml.in \
  SPAMO_config.xml.in \
  TOMTOM_config.xml.in

BUILT_CONFIG_FILES = \
  AME_$(VERSION).xml \
  CENTRIMO_$(VERSION).xml \
  DREME_$(VERSION).xml \
  FIMO_$(VERSION).xml \
  GLAM2_$(VERSION).xml \
  GLAM2SCAN_$(VERSION).xml \
  GOMO_$(VERSION).xml \
  MAST_$(VERSION).xml \
  MCAST_$(VERSION).xml \
  MEMECHIP_$(VERSION).xml \
  MEME_$(VERSION).xml \
  MOMO_$(VERSION).xml \
  SPAMO_$(VERSION).xml \
  TOMTOM_$(VERSION).xml

SOURCE_ETC_FILES = \
  constants.js.in \
  constants.xsl.in \
  logging.conf.in \
  xsl-wrapper.in

NORMAL_ETC_FILES = \
  alphabet_utilities.js \
  ame_icon.png \
  ame_template.html \
  ame_template.css \
  ame_template.js \
  beadstring.xsl \
  block_diagram.css \
  block_diagram.js \
  centrimo_graph.js \
  centrimo_icon.png \
  centrimo_template.html \
  centrimo_template.css \
  centrimo_template.js \
  cisml.css \
  cisml-to-gff.xsl \
  cisml-to-gff3.xsl \
  cisml-to-html.xsl \
  cisml-to-text.xsl \
  delay_draw.js \
  dreme_icon.png \
  dreme_template.html \
  dreme_template.css \
  dreme_template.js \
  eps_context.js \
  fimo_icon.png \
  fimo-to-gff3.xsl \
  fimo-to-html.xsl \
  fimo-to-text.xsl \
  generate_js_wrapped_template.pl \
  generate_xsl_icons.pl \
  gomo-to-html.xsl \
  hmm-text.xsl \
  job_status.tmpl \
  mast_icon.png \
  mast_template.css \
  mast_template.html \
  mast_template.js \
  mcast_icon.png \
  mcast_template.css \
  mcast_template.html \
  mcast_template.js \
  meme.css \
  meme_icon.png \
  meme_template.html \
  meme_template.css \
  meme_template.js \
  memechip_icon.png \
  meme-chip_template.html \
  meme-chip_template.css \
  meme-chip_template.js \
  motif_logo.js \
  new_icon.png \
  prior1.plib \
  prior30.plib \
  spamo_graph.js \
  spamo_icon.png \
  spamo_template.html \
  spamo_template.css \
  spamo_template.js \
  template.eps \
  tomtom_icon.png \
  tomtom_template.html \
  tomtom_template.css \
  tomtom_template.js \
  utilities.js

BUILT_ETC_FILES = \
  constants.js \
  constants.xsl \
  icons.xsl \
  logging.conf \
  meme.css.xsl \
  motif_logo_template.js \
  utilities.js.xsl

%.xsl: %
	$(SED) 's/tname/$?/' xsl-wrapper.in | $(SED) '5r $?' > $@

%_$(VERSION).xml: %_config.xml.in
	$(SED) -e 's%@SERVICE_VERSION@%$(VERSION)%' -e 's%@SERVICE_DIR@%$(bindir)%' $? > $@

SED_URLS = -e 's%@SITE_URL@%$(URL)%' -e 's%@AMIGO_URL@%$(GOSERVER_URL)%'

constants.js: constants.js.in
	$(SED) $(SED_URLS) $? > $@

constants.xsl: constants.xsl.in
	$(SED) $(SED_URLS) $? > $@

icons.xsl: ../doc/images/gomo_icon.png
	$(PERL) generate_xsl_icons.pl $? > $@

motif_logo_template.js: template.eps
	$(PERL) generate_js_wrapped_template.pl motif_logo_template $? > $@

logging.conf: logging.conf.in 
	$(SED) $(SED_URLS) $? > $@

etcdir=$(sysconfdir)

etc_DATA = $(DATABASE_FILES) $(NORMAL_ETC_FILES) $(BUILT_ETC_FILES)

if WEBSERVICE
deploydir=$(OPAL_DEPLOY)

deploy_DATA = $(BUILT_CONFIG_FILES)
endif
 
CLEANFILES = $(BUILT_CONFIG_FILES) $(BUILT_ETC_FILES)

EXTRA_DIST = $(DATABASE_FILES) $(SOURCE_ETC_FILES) $(NORMAL_ETC_FILES) $(SOURCE_CONFIG_FILES)


