##
## Process this file with automake to produce Makefile.in
##

#
# To MAKE in DEBUG mode (it will write out lots of status info):
#	make CFLAGS=" -DDEBUG_PARALLEL"
#
# To MAKE in DEBUG mode and force the serial compile so GDB will work:
#	make MPICC=gcc MPIFLAGS="" MPI_CMD="" CFLAGS=" -ggdb -O0 -DDEBUG"

# Note: All programs that can be compiled in parallel mode are
# compiled in this directory.  We only make the parallel versions
# if we can.  Otherwise, we make the serial versions.

# Note: If we are only making the serial versions, MPICC is set to 
# the serial compiler and MPIFLAGS is empty.
CC = $(MPICC)

bin_PROGRAMS = meme momo
meme_LDADD = ../libcommon.la $(LIBXML2_LIBS)
meme_CFLAGS = -I.. $(MPIFLAGS) $(LIBXML2_CFLAGS) -DMPI_CMD="${MPI_CMD}" -DMEME_MAIN

meme_SOURCES = \
  ../banner.c \
  ../branching_search.c \
  ../calculate_p_y.c \
  ../discretize.c \
  ../display.c \
  ../dpalign.c \
  ../em.c  \
  ../exec_parallel.c \
  ../fisher_exact.c \
  ../init.c  \
  ../justlike.c \
  ../likelihood.c \
  ../llr.c \
  ../logs.c \
  ../meme.c \
  ../meme_util.c \
  ../message.c \
  ../mp.c \
  ../oops.c \
  ../partition.c \
  ../prior.c  \
  ../psp.c \
  ../read_seq_file.c \
  ../seed.c \
  ../seed_diffs.c \
  ../seq2theta.c \
  ../sp_matrix.c \
  ../starts.c \
  ../subseq7.c \
  ../tcm.c \
  ../ushuffle.c \
  ../verbosity.c

momo_CFLAGS = $(AM_CFLAGS) $(LIBXML2_CFLAGS)
momo_LDADD = ../libcommon.la $(LIBXML2_LIBS)
momo_SOURCES = ../momo.c \
  ../momo-output.c \
  ../momo-algorithm.c \
  ../momo-simple.c \
  ../momo-motifx.c \
  ../momo-modl.c \
  ../momo-input.c \
  ../banner.c \
  ../branching_search.c \
  ../calculate_p_y.c \
  ../discretize.c \
  ../display.c \
  ../dpalign.c \
  ../em.c  \
  ../exec_parallel.c \
  ../fisher_exact.c \
  ../init.c  \
  ../justlike.c \
  ../likelihood.c \
  ../llr.c \
  ../logs.c \
  ../meme.c \
  ../meme_util.c \
  ../message.c \
  ../mp.c \
  ../oops.c \
  ../partition.c \
  ../prior.c  \
  ../psp.c \
  ../read_seq_file.c \
  ../seed.c \
  ../seed_diffs.c \
  ../seq2theta.c \
  ../sp_matrix.c \
  ../starts.c \
  ../subseq7.c \
  ../tcm.c \
  ../ushuffle.c \
  ../verbosity.c

check:
	cp -f meme ../meme; \
	cp -f momo  ../momo;
