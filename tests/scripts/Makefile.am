##
## Process this file with automake to produce Makefile.in
##

AUTOMAKE_OPTIONS = no-dependencies

test_driver: test_driver.pl.in
	$(SED) -e 's%@WHICHPERL@%$(PERL)%' $< > $@
	chmod +x $@

TEST_FILES = ama.test ame.test centrimo.test \
  create_priors.test dreme.test fasta-center.test \
  fasta-dinucleotide-shuffle.test fimo.test glam2.test \
  glam2scan.test gomo.test mast.test \
  mcast.test meme.test meme-chip.test momo.test motif-in.test \
	motiph.test pspgen.test qvalue.test spamo.test tomtom.test

# tests for programs that are not built by default
OBSOLETE_TEST_FILES = clustalw2fasta.test fasta-io.test mhmm_core.test \
  mhmm-io.test mhmms.test mhmmscan.test shadow.test 

# Broken tests
XFAIL_TESTS = 

check_SCRIPTS = test_driver

TEST_LOG_DRIVER = $(PERL) test_driver

TESTS = $(TEST_FILES)

EXTRA_DIST = $(DATA_FILES) $(TEST_FILES) $(OBSOLETE_TEST_FILES) test_driver.pl.in

CLEANFILES = test_driver
