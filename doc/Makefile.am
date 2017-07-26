#####################################################################
## FILE: Makefile.am
# AUTHOR: William Stafford Noble, Charles E. Grant
# CREATE DATE: 2003
# PROJECT: Meta-MEME
# COPYRIGHT: 2003, UW
#####################################################################

## Automake uses this file to generate a Makefile

#####################################################################
#  This is the Makefile for the documentation directory.
#####################################################################
SUBDIRS = examples images css js

copyright.html: copyright.html.in ../COPYING Makefile
	$(SED) -e "/INSERT_COPYING_FILE/r ../COPYING" -e "/INSERT_COPYING_FILE/d" copyright.html.in > copyright.html


# Files Included in distribution and installed.
NORMAL_FILES = \
  ../etc/prior30.plib \
  alphabets.html \
  alphabet-format.html \
  alphtype.html \
  ama.html \
  ama-qvalues.html \
  ame.html \
  ame-tutorial.html \
  authors.html \
  beadstring.html \
  beeml2meme.html \
  bfile-format.html \
  centrimo.html \
  centrimo-tutorial.html \
  ceqlogo.html \
  chen2meme.html \
  cite.html \
  clustalw2fasta.html \
  clustalw2phylip.html \
  clustalw-format.html \
  compute-prior-dist.html \
  compute-uniform-priors.html \
  copyright.html \
  create-priors.html \
  dmix-format.html \
  download.html \
  draw-mhmm.html \
  dreme.html \
  dreme-tutorial.html \
  dust.html \
  elm2meme.html \
  fasta-fetch.html \
  fasta-format.html \
  fasta-get-markov.html \
  fasta-grep.html \
  fasta-hamming-enrich.html \
  fasta-io.html \
  fasta-make-index.html \
  fasta-most.html \
  fasta-shuffle-letters.html \
  fasta-subsample.html \
  fasta-unique-names.html \
  fimo.html \
	fimo-tutorial.html \
  fisher_exact.html \
  fitevd.html \
  gendb.html \
  general-faq.html \
  getsize.html \
  glam2.html \
  glam2format.html \
  glam2mask.html \
  glam2scan.html \
  glam2_alphabet.html \
  glam2_tut.html \
  GLAM2_method.pdf \
  godag-format.html \
  gomo.html \
  gomo-format.html \
  gomo_highlight.html \
  install.html \
  iupac.html \
  iupac2meme.html \
  jaspar2meme.html \
  log-hmm.html \
  mast.html \
  matrix2meme.html \
  mcast.html \
  meme.html \
  meme2alph.html \
  meme-chip.html \
  meme-format.html \
  meme-get-motif.html \
  meme-rename.html \
  meme2images.html \
  meme2meme.html \
  mhmme.html \
  mhmm.html \
  mhmm2html.html \
  mhmm-format.html \
  mhmmscan.html \
  mhmmscan-format.html \
  mhmms.html \
  momo.html \
  motif-shuffle-columns.html \
  motif-consensus.html \
  motif_conversion.html \
  motiph.html \
  nmica2meme.html \
  overview.html \
  pmp_bf.html \
  priority2meme.html \
	psm-format.html \
  psp-format.html \
  psp-gen.html \
  purge.html \
  qvalue.html \
  reconcile-tree-alignment.html \
  reduce-alignment.html \
  release-notes.html \
  remove-alignment-gaps.html \
  rna2meme.html \
  scpd2meme.html \
  scripting.html \
  shadow.html \
  sites2meme.html \
  spamo.html \
  taipale2meme.html \
  tamo2meme.html \
  tomtom.html \
  transfac2meme.html \
  transfac-format.html \
  uniprobe2meme.html \
  update-sequence-db.html \
  ismb94.pdf

memedocdir = $(prefix)/share/doc
memedoc_DATA =  $(NORMAL_FILES)

# force the generation of copyright.html when making a website
BUILT_SOURCES = copyright.html

EXTRA_DIST = copyright.html.in $(NORMAL_FILES)
