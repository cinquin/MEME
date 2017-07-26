##
## Process this file with automake to produce Makefile.in
##

##############################################
# some examples are also used in tests: check
# scripts/metatest.pl.in before changing

AUTOMAKE_OPTIONS = no-dependencies


SUBDIRS = \
  compute_prior_dist_example_output_files \
  sample_opal_scripts

EXAMPLES = \
  sample-dna-aln.glam2 \
  sample-protein-aln.glam2 \
  sample-dna-motif.meme \
  sample-protein-motif.meme \
  small-tide-mod.txt

# load the EXAMPLE_OUTPUT_FILES variable
include examples.mk

if ! WEBSITE
exdir = $(prefix)/share/doc/examples
nobase_dist_ex_DATA = $(EXAMPLES) $(EXAMPLE_OUTPUT_FILES)
endif

EXTRA_DIST = examples.mk $(EXAMPLES) $(EXAMPLE_OUTPUT_FILES)
