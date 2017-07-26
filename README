***Introduction***           

The  MEME suite provides online tools for discovering and using protein and 
DNA sequence motifs.  A motif is a pattern of nucleotides or amino acids that 
appears repeatedly in a group of related  DNA or protein sequences.  The MEME 
suite represents motifs as position-dependent scoring matrices. 

The MEME suite consists of programs which allow you to: 

* meme      - for discovery of motifs shared by a group of sequences, 
* mast      - for search of databases for sequences containing these motifs,
* tomtom    - for searching databases of motifs for similar motifs,
* gomo      - for finding Gene Ontology terms linked to the motifs,
* glam2     - for discovery of gapped motifs,
* glam2scan - for scanning sequences with gapped motifs,
* fimo      - for scanning sequences with motifs,
* mcast     - for finding motif clusters,
* meme-chip - for analysis of large DNA datasets like ChIPseq output,
* spamo     - for finding motif complexes by analysing motif spacing,
* dreme     - for discovery of short regular expression motifs,

and that's just the web enabled tools.

You can download the C source code for MEME suite from
http://meme-suite.org/doc/download.html

You can also use the MEME suite via its website at http://meme-suite.org .

***Citing***

To cite the full MEME suite, please cite:  
> Timothy L. Bailey, Mikael Bodén, Fabian A. Buske, Martin Frith,   
  Charles E. Grant, Luca Clementi, Jingyuan Ren, Wilfred W. Li,   
  William S. Noble, "MEME SUITE: tools for motif discovery and searching",   
  Nucleic Acids Research, 37:W202-W208, 2009.  

To cite individual tools, please check the citation page:  
    http://meme-suite.org/doc/cite.html  


***Installation***       

See `doc/install.html` for operating system requirements, prerequisite
software, and installation instructions.

***Documentation***

Documentation is available online at http://meme-suite.org/doc/overview.html
otherwise look in the `doc/` folder for `overview.html` as a place to start.  
If you did installation with the `--enable-web` switch, the html documentation
will be installed with the website at `<website>/meme_<version>/doc/overview.html`. 

***Problems and comments*** 

Please address any problems or comments to:  
    meme-suite@uw.edu 
or  
    https://groups.google.com/forum/#!forum/meme-suite  

***Release Notes***

See file `<distribution-path>/doc/release-notes.html`  
or after a basic install see `<install-path>/doc/release_notes.html`  
or for a website install see `<website>/meme_<version>/doc/release-notes.html`  

***Developers Notes***

To prepare a new release.

1. Clone from bitbucket (you have to be granted access to this mercurial repository):  
> hg clone ssh://hg@bitbucket.org/tlbailey/meme

2. Create a release branch (Note: the version number is set in configure.ac using variable AM_INIT_AUTOMAKE.):
> hg branch meme_VERSION  
  hg ci -m "Create release branch for VERSION"

3. Once you are ready to create a release candidate tag the current revision
(this is used to determine the release date and release revision):  
> hg tag meme_VERSION_0  

    3.1. If you later have to create another release candidate after applying some
    patches then move the tag so that the release date and release revision are correct:  
    > hg tag -f meme_VERSION_0

4. Create the example output files by building a copy of the MEME Suite, ensuring
it is first in your path, and running create_examples.pl:
> cd meme/doc/examples  
  ./create_examples.pl  
  scp examples.tgz http://meme-suite.org/meme-software/example-output/examples.tgz # ok so this line doesn't really work but you should get the intent

5. To create a distribution tar ball meme_VERSION.tar.gz.
> cd meme  
  hg purge --all  
  ./bootstrap  
  ./configure --enable-web --enable-opt --enable-build-libxml2 --enable-build-libxslt  
  make dist
