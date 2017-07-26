
MemeMenu.blurbs = {
  "memesuite": "The MEME Suite allows you to discover novel motifs in collections of unaligned nucleotide or protein sequences, and to perform a wide variety of other motif-based analyses.  </br></br> It provides <em>motif discovery algorithms</em> using both probabilistic (MEME) and discrete models (DREME), which have complementary strengths.  It also allows discovery of motifs with <em>arbitrary insertions and deletions</em> (GLAM2).  MEME-ChIP performs a <em>comprehensive motif-based analysis of ChIP-seq</em> and other large sequence datasets.  </br></br>In addition to motif discovery, you can also perform several kinds of motif enrichment analysis using compendia of known motifs:  measuring the <em>positional enrichment</em> of sequences for known motifs (CentriMo), measuring the <em>relative enrichment</em> between two sets of sequences (AME), detecting <em>preferred spacings</em> between pairs of motifs (SpaMo) and identifying the <em>biological roles</em> of motifs (GOMo).  </br></br>The MEME Suite also provides tools for <em>scanning sequences</em> for matches to motifs (FIMO, MAST and GLAM2Scan), and for scanning for clusters of motifs (MCAST).  </br></br>Finally, you can <em>compare motifs</em> to known motifs using Tomtom.</br></br>Position the cursor on the flow chart below and/or click for more information on the tools in the MEME suite and information on their inputs and outputs.  Clicking on any of the tool icons also takes you to their input forms.  </br></br><b>Click</b> here for more information on the MEME Suite.",

  "input_sequences": "You can use the MEME Suite tools to discover novel (<em>Motif Discovery</em>) or known (<em>Motif Enrichmment</em>) sequence motifs in sets of related DNA, RNA or protein sequences.  </br></br>You can also input sets of sequences and scan them for occurrences of motifs (<em>Motif Scanning</em>).  </br></br>Details on the format of your sequences are given under <em>FASTA Sequence</em> in the <em>File Format Reference</em> menu on the left, or just by <b>clicking</b> here.",

  "input_motifs": "You can input your own motifs to MEME Suite tools to see if they are enriched in your sequences (<dfn>Motif Enrichment</dfn>), to find out where they occur in known sequences (<dfn>Motif Scanning</dfn>), or to see if they are similar to known motifs (<dfn>Motif Comparison</dfn>).  </br></br>Many motif formats are supported including count matrix, postition weight matrix, aligned sites, and consensus sequence.  </br></br>You can also directly input the motifs contained in the output of the MEME Suite Motif Discovery tools, or a simplified version of that format, which is detailed under <dfn>MEME Motif</dfn> in the <dfn>File Format Reference</dfn> menu on the left, or just by <b>clicking</b> here.", 

  "motif_databases": "The MEME Suite provides a large number of databases of known motifs that you can use with the <dfn>Motif Enrichment</dfn> and <dfn>Motif Comparison</dfn> tools. </br></br>The motif databases are also available for you to download and use on your own computer under <dfn>Download MEME Suite and Databases</dfn> in the <dfn>Download & Install</dfn> menu on the left. </br></br><b>Click</b> here to see descriptions of the available motif databases.",

  "go_databases": "The MEME Suite provides databases for several organisms for use with the GOMo motif enrichment tool. </br></br>These databases consist of the promoter sequences for all annotated genes in the organism, plus a table listing the Genome Ontology (GO) terms associated with each gene. </br></br><b>Click</b> here to see descriptions of the available GOMo sequence databases. </br></br>The GOMo databases are also available for you to download and use on your own computer under <dfn>Download MEME Suite and Databases</dfn> in the <dfn>Download & Install</dfn> menu on the left.",

  "discovered_motifs": "The MEME Suite <dfn>Motif Discovery</dfn> tools output their results as interactive HTML files. </br></br>You can see an example output from each of the tools under <dfn>Motif Discovery</dfn> in the <dfn>Sample Outputs</dfn> menu on the left.",

  "enriched_motifs": "The MEME Suite <dfn>Motif Enrichment</dfn> tools output their results as interactive HTML files. </br></br>You can see an example output from each of the tools under <dfn>Motif Enrichment</dfn> in the <dfn>Sample Outputs</dfn> menu on the left.",

  "annotated_motifs": "The MEME Suite <dfn>GOMo</dfn> tool outputs its results as an interactive HTML file. </br></br>You can see an example <dfn>GOMo</dfn> output under <dfn>Motif Enrichment</dfn> in the <dfn>Sample Outputs</dfn> menu on the left, or by <b>clicking</b> here.",

  "sequence_databases": "The MEME Suite provides a large number of protein and DNA sequence databases.  These include the proteomes, genomes and sets of promoters for many organisms.</br></br><b>Click</b> here to see descriptions of the available sequence databases.",
  "motif_discovery": "<i>De novo</i> discovery of motifs in set(s) of sequences.",
  "motif_enrichment": "Analyze sequence set(s) for enrichment of known motifs or motifs you provide.",
  "motif_scanning": "Find matches to motif(s) in sequences.",
  "motif_comparison": "Compare query motif(s) to known motifs.",

  "momo": "<dfn>MoMo</dfn> identifies sequence motifs associated with different types of protein post-translational modifications (PTMs). The program takes as input a collection of PTMs identified using protein mass spectrometry.  For each distinct type of PTM, the momo command uses one of three algorithms to constructs a motifs representing amino acid preferences flanking the modification site.</footer>",
  
  "meme": "<dfn>MEME</dfn> discovers novel, <em>ungapped</em> motifs (recurring, fixed-length patterns) in your sequences<span class=\"sample_output\"> (<a href=\"../doc/examples/meme_example_output_files/meme.html\">sample output</a> from <a href='http://meme-suite.org/meme-software/example-datasets/crp0.fna'>sequences</a>)</span>. <dfn>MEME</dfn> splits variable-length patterns into two or more separate motifs. <span class=\"manual_link\">See this <a href=\"../doc/meme.html?man_type=web\">Manual</a> for more information.</span> <footer><br><br>A motif is a sequence pattern that occurs repeatedly in a group of related sequences. MEME represents motifs as position-dependent letter-probability matrices which describe the probability of each possible letter at each position in the pattern.  Individual MEME motifs do not contain gaps. Patterns with variable-length gaps are split by MEME into two or more separate motifs. <br><br>MEME takes as input a group of sequences and outputs as many motifs as requested. MEME uses statistical modeling techniques to automatically choose the best width, number of occurrences, and description for each motif. <br><br>MEME on the web can take a second (control) set of input sequences and then discovers motifs that are enriched in the primary set relative to the control set.  This is called discriminative motif discovery.</footer>",

  "dreme": "<dfn>DREME</dfn> discovers <em>short</em>, <em>ungapped</em> motifs (recurring, fixed-length patterns) that are <em>relatively</em> enriched in your sequences compared with shuffled sequences or your control sequences<span class=\"sample_output\"> (<a href=\"../doc/examples/dreme_example_output_files/dreme.html\">sample output</a> from <a href='http://meme-suite.org/meme-software/example-datasets/Klf1.fna'>sequences</a>)</span>.<span class=\"manual_link\"> See this <a href=\"../doc/dreme.html?man_type=web\">Manual</a> or this <a href=\"../doc/dreme-tutorial.html\">Tutorial</a> for more information.</span>  <footer><br><br>DREME (Discriminative Regular Expression Motif Elicitation) finds relatively short motifs (up to 8 positions) fast. The input to DREME is one or two sets of sequences.  The control sequences should be approximately the same length as the primary sequences. If you do not provide a control set, the program shuffles the primary set to provide a control set. The program uses Fisher's Exact Test to determine significance of each motif found in the postive set as compared with its representation in the control set, using a significance threshold that may be set on the command line.<br><br>DREME achieves its high speed by restricting its search to regular expressions based on the symbols available in the alphabet, and by using a heuristic estimate of generalised motifs' statistical significance.</footer>",

  "memechip": "<dfn>MEME-ChIP</dfn> performs <em>comprehensive motif analysis</em> (including motif discovery) on LARGE (50MB maximum) sets of sequences (typically <em>nucleotide</em>) such as those identified by ChIP-seq or CLIP-seq experiments<span class=\"sample_output\"> (<a href=\"../doc/examples/memechip_example_output_files/meme-chip.html\">sample output</a> from <a href='http://meme-suite.org/meme-software/example-datasets/Klf1.fna'>sequences</a>)</span>.<span class=\"manual_link\"> See this <a href=\"../doc/meme-chip.html?man_type=web\">Manual</a> for more information.</span>  <footer><br><br>MEME-ChIP can: <ol> <li>discover novel DNA-binding motifs (with MEME and DREME),</li> <li>determine which motifs are most centrally enriched (with CentriMo),</li> <li>analyze them for similarity to known binding motifs (with Tomtom), and</li> <li>automatically group significant motifs by similarity,</li><li>perform a motif spacing analysis (with SpaMo), and,</li><li>create a GFF file for viewing each motif's predicted sites in a genome browser.</li></ol>It is worth noting that MEME-ChIP is <b>not</b> a motif scanner, but the motifs discovered by it can be used by FIMO or MAST to scan for motif sites.</footer>",

  "glam2": "<dfn>GLAM2</dfn> discovers novel, <em>gapped</em> motifs (recurring, variable-length patterns) in your <em>DNA</em> or <em>protein</em> sequences<span class=\"sample_output\"> (<a href=\"../doc/examples/glam2_example_output_files/glam2.html\">sample output</a> from <a href='http://meme-suite.org/meme-software/example-datasets/At.faa'>sequences</a>)</span>.<span class=\"manual_link\"> See this <a href=\"../doc/glam2.html?man_type=web\">Manual</a> for more information.</span>  <footer><br><br>GLAM2 is a program for finding motifs in sequences, typically amino-acid or nucleotide sequences. The main innovation of GLAM2 is that it allows insertions and deletions in motifs.  <br><br>GLAM2 is a child of <a href=\"http://zlab.bu.edu/glam/\">GLAM</a> and a sibling of <a href=\"http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html.ncbi/index/software.html#7\">A-GLAM</a>.  'GLAM' previously to stood for 'Gapless Local Alignment of Multiple sequences', but it now stands for 'Gapped Local Alignment of Motifs'.  <br><br>GLAM2 uses a <a href=\"http://en.wikipedia.org/wiki/Simulated_annealing\">simulated annealing</a> algorithm, with a temperature parameter. At high temperatures, GLAM2 only slightly favours changes that increase the alignment's score, and at low temperatures it strongly favours such changes. Thus, at high temperatures the score will be optimized too slowly, but at low temperatures the algorithm will get frozen in a local optimum. The strategy, then, is to start with a high temperature and reduce it as slowly as possible.  <br><br>Refer to the <a href=\"glam2_tut.html\">GLAM2 tutorial</a>, which uses the command-line version of GLAM2, if you need more information.</footer>",

  "ame": "<dfn>AME</dfn> identifies <em>known</em> or <em>user-provided</em> motifs that are <em>relatively</em> enriched in your sequences compared with shuffled versions of your sequences (web version only) or with your control sequences<span class=\"sample_output\"> (<a href=\"../doc/examples/ame_example_output_files/ame.html\">sample output</a> from <a href='http://meme-suite.org/meme-software/example-datasets/crp0.fna'>sequences</a>, <a href='http://meme-suite.org/meme-software/example-datasets/lex0.fna'>control sequenes</a> and <a href='http://meme-suite.org/meme-software/example-datasets/dpinteract.meme'>motifs</a>)</span>. <dfn>AME</dfn> treats motif occurrences the same, regardless of their locations within the sequences. (The command-line version of AME can also look for enrichment in a single file where the sequences are sorted by some biological signal.) <span class=\"manual_link\">See this <a href=\"../doc/ame.html?man_type=web\">Manual</a> or this <a href=\"../doc/ame-tutorial.html\">Tutorial</a> for more information.</span> <footer><br><br>AME (Analysis of Motif Enrichment) scores a set of sequences given a set of binding motifs, treating each subsequence (and its reverse complement for complementable alphabets) in the sequence as a possible match to the motif. The web version of AME provides a large number of motif databases, or you can upload your own motifs in <a href=meme-format.html>MEME Motif Format</a>, or enter them by typing in many additional formats.</footer>", 
  "centrimo": "<dfn>CentriMo</dfn> identifies <em>known</em> or <em>user-provided</em> motifs that show a significant preference for particular <em>locations</em> in your sequences<span class=\"sample_output\"> (<a href=\"../doc/examples/centrimo_example_output_files/centrimo.html\">sample output</a> from <a href='http://meme-suite.org/meme-software/example-datasets/mm9_tss_500bp_sampled_1000.fna'>sequences</a> and <a href='http://meme-suite.org/meme-software/example-datasets/some_vertebrates.meme'>motifs</a>)</span>. <dfn>CentriMo</dfn> can also show if the local enrichment is significant <em>relative</em> to control sequences. <span class=\"manual_link\">See this <a href=\"../doc/centrimo.html?man_type=web\">Manual</a> or this <a href=\"../doc/centrimo-tutorial.html\">Tutorial</a> for more information.</span> <footer><br><br>CentriMo takes a set of motifs and a set of equal-length sequences and plots the positional distribution of the best match of each motif.  The motifs are typically compendia of DNA- or RNA-binding motifs, and the sequences might be: 500 bp sequences aligned on ChIP-seq peaks or summits; 300 bp sequences centered on sets of transcription start sites or translation start sites; sequences aligned on splice-junctions; etc.  The web version of CentriMo provides a large number of motif databases, or you can upload your own motifs in <a href=meme-format.html>MEME Motif Format</a>, or enter them by typing in many additional formats. <br><br>CentriMo also computes the 'local enrichment' of each motif by counting the number of times its best match in each sequence occurs in a local region and applying a statistical test to see if the local enrichment is signficant.  By default, CentriMo examines only regions <b>centered</b> in the input sequences, but CentriMo will compute the enrichment of <b>all</b> regions if requested.  CentriMo uses the binomial test to compute the significance of the number of sequences where the best match occurs in a given region, assuming a uniform prior over best match positions.  CentriMo reports the location and significance of the best region for each motif.  <br><br> CentriMo can also perform comparative enrichment, reporting the relative enrichment of the best region in a second, control set of sequences if you provide a set of control sequences.  CentriMo choses regions based on their significance in the primary set of sequences, and then it uses Fisher's exact test to evaluate the significance of the number of best matches in the region in the primary set compared with the number of best matches in the same region in the control set of sequences.</footer>",

  "spamo": "<dfn>SpaMo</dfn> identifies significantly enriched <em>spacings</em> in a set of sequences between a <em>primary motif</em> and each motif in a set of <em>secondary motifs</em><span class=\"sample_output\"> (<a href=\"../doc/examples/spamo_example_output_files/spamo.html\">sample output</a> from <a href='http://meme-suite.org/meme-software/example-datasets/Klf1.fna'>sequences</a> and <a href='http://meme-suite.org/meme-software/example-datasets/Klf1.meme'>motif</a>)</span>. Typically, the input sequences are centered on ChIP-seq peaks, and are each 500bp long.  <span class=\"manual_link\">See this <a href=\"../doc/spamo.html?man_type=web\">Manual</a> for more information.</span>  <footer><br><br>The name SpaMo stands for 'Spaced Motif' analysis.  Its inputs are a set of many short sequences, a primary motif and one or more databases of secondary motifs. It searches for the strongest primary motif binding site and then searches in the area around it looking for the strongest secondary motif binding site. The relative spacings of the primary and secondary motif in all the sequences is tallied and the probability of the close spacings happening by chance is calculated.  <br><br> After all the calculations are done SpaMo, outputs the non-redundant secondary motifs in order of significance provided they had a bin that passed the significance threshold. Similar secondary motifs are grouped together and listed in order of significance on the secondary motif they were redundant to. If the bin size is 1 (as it is with the web version of SpaMo), then an alignment for each of the similar secondary motifs is created.  <br><br>The web version of SpaMo provides a large number of secondary motif databases, or you can upload your own motifs in <a href=meme-format.html>MEME Motif Format</a>, or enter them by typing in many additional formats. <br><br>It is important to note that the locations searched for the primary and secondary motifs are dependant on a SpaMO parameter called the \"margin\", which is 150bp by default. The secondary motif must be contained within an area limited by the margin size on either side of the strongest primary motif site. Additionally the primary motif cannot be found closer than the margin size to the edge of the input sequence. <br><br><figure> <figcaption>Restrictions on primary motif locations.</figcaption> <canvas id=\"primary_site_diagram\"></canvas> </figure> <figure> <figcaption>Restrictions on secondary motif locations.</figcaption> <canvas id=\"secondary_site_diagram\" onclick=\"draw_diagram2()\"></canvas> </figure></footer>",

  "gomo": "<dfn>GOMo</dfn> scans all <em>promoters</em> using <em>nucleotide motifs</em> you provide to determine if any motif is significantly associated with genes linked to one or more Genome Ontology (GO) terms<span class=\"sample_output\"> (<a href=\"../doc/examples/gomo_example_output_files/gomo.html\">sample output</a> from <a href='http://meme-suite.org/meme-software/example-datasets/dpinteract_subset.meme'>motifs</a> and E.coli K12 database)</span>. The significant GO terms can suggest the <em>biological roles</em> of the motifs.  <span class=\"manual_link\">See this <a href=\"../doc/gomo.html?man_type=web\">Manual</a> for more information.</span>  <footer><br><br>The name GOMo stands for 'Gene Ontology for Motifs'.  The program searches in a set of ranked genes for enriched GO terms associated with high-ranking genes. The genes can be ranked, for example, by applying a motif scoring algorithms on their upstream sequence.  <br><br>GOMo computes a score for each GO term by estimating the Mann-Whitney rank-sum <i>p</i>-value of the GO term's genes, and can combine the <i>p</i>-values for multiple species. It then estimates final <i>p</i>-values for each GO-term empirically by shuffling the gene identifiers in the ranking (ensuring consistency across species) to generate scores from the null hypothesis. Then <i>q</i>-values are derived from these <i>p</i>-values following the method of Benjamini and Hochberg (where '<i>q</i>-value' is defined as the minimal false discovery rate at which a given GO-term is deemed significant).<br><br>The program reports all GO terms that receive <i>q</i>-values smaller than a specified threshold, outputting a GOMo score with emprically calculated <i>p</i>-values and <i>q</i>-values for each.</footer>",

  "fimo": "<dfn>FIMO</dfn> scans a sequence database for <em>individual matches</em> to each of the motifs you provide<span class=\"sample_output\"> (<a href=\"../doc/examples/fimo_example_output_files/fimo.html\">sample output</a> for <a href='http://meme-suite.org/meme-software/example-datasets/some_vertebrates.meme'>motifs</a> and <a href='http://meme-suite.org/meme-software/example-datasets/mm9_tss_500bp_sampled_1000.fna'>sequences</a>)</span>.<span class=\"manual_link\"> See this <a href=\"../doc/fimo.html?man_type=web\">Manual</a> or this <a href=\"../doc/fimo-tutorial.html\">Tutorial</a> for more information.</span>  <footer><br><br>The name FIMO stands for 'Find Individual Motif Occurrences'.  The program searches a database of sequences for occurrences of known motifs, treating each motif independently.  Motifs must be in <a href=meme-format.html>MEME Motif Format</a>. The web version of FIMO also allows you to type in motifs in additional formats. <br><br> You can define the statistical threshold (<i>p</i>-value) for motifs and whether FIMO scans just the given sequences or their reverse complements (where applicable), too.</footer>", 

  "ptm": "<dfn>PTM Portal</dfn> discovers motifs from inputs of peptide-spectrum matches.",

  "mast": "<dfn>MAST</dfn> searches sequences for matches to a set of motifs, and sorts the sequences by the <em>best combined match</em> to all motifs<span class=\"sample_output\"> (<a href=\"../doc/examples/mast_example_output_files/mast.html\">sample output</a> for <a href='http://meme-suite.org/meme-software/example-datasets/adh.meme'>motifs</a> and <a href='http://meme-suite.org/meme-software/example-datasets/adh.faa'>sequences</a>).</span><span class='manual_link'> See this <a href=\"../doc/mast.html?man_type=web\">Manual</a> for more information.</span>  <footer><br><br>The name MAST stands for 'Motif Alignment and Search Tool'.  A motif is a sequence pattern that occurs repeatedly in a group of related sequences. Motifs are represented as position-dependent probability (or scoring) matrices that describe the probability (or score) of each possible letter at each position in the pattern. Individual motifs may not contain gaps. Patterns with variable-length gaps must be split into two or more separate motifs before being submitted as input to MAST.  <br><br>MAST takes as input a file containing the descriptions of one or more motifs and searches a sequence database that you select for sequences that match the motifs.  Motifs must be in <a href=meme-format.html>MEME Motif Format</a>. The web version of MAST also allows you to type in motifs in additional formats.  <br><br>MAST can search sequences with motifs of the same alphabet but as a special feature it can also search protein sequences with DNA motifs.  You can control how (or if) MAST scans the reverse complement strand (if applicable).  You can also set thresholds for overall (combined match) significance and individual motif site significance.</footer>",

  "mcast": "<dfn>MCAST</dfn> searches sequences for <em>clusters of matches</em> to one or more <em>nucleotide</em> motifs<span class=\"sample_output\"> (<a href=\"../doc/examples/mcast_example_output_files/mcast.html\">sample output</a> for <a href='http://meme-suite.org/meme-software/example-datasets/crp0.meme'>motifs</a> and <a href='http://meme-suite.org/meme-software/example-datasets/crp0.fna'>sequences</a>)</span>.<span class='manual_link'> See this <a href=\"../doc/mcast.html?man_type=web\">Manual</a> for more information.</span>  <footer><br><br>The name MCAST stands for 'Motif Cluster Alignment Search Tool'.  MCAST searches a sequence database for statistically significant clusters of <b>non-overlapping</b> occurrences of a given set of motifs.  It is primarily designed for search genomic DNA databases. Motifs must be in <a href=meme-format.html>MEME Motif Format</a>. The web version of MCAST also allows you to type in motifs in additional formats.  <br><br>A motif 'hit' is a sequence position that is sufficiently similar to a motif in the query, where the score for a motif at a particular sequence position is computed without gaps. To compute the <i>p-</i>value of a motif score, MCAST assumes that the sequences in the database were generated by a 0-order Markov process.  To be considered a hit, the <i>p</i>-value of the motif alignment score must be less than the significance threshold, which you can specify, Note that MCAST searches for hits on both strands of the sequences.  <br><br>A cluster of non-overlapping hits is called a 'match'. You can specify the maximum allowed distance between the hits in a match.  Two hits separated by more than the maximum allowed gap will be reported in separate matches.  <p style=\"margin-bottom:0\">The <i>p-</i>value of a hit is converted to a 'p-score' in order to compute the total score of the match it participates in. The p-score for a hit with <i>p-</i>value <span class=\"pdat\">p</span> is</p> <div class=\"indent\"> S = -log<sub>2</sub>(<span class=\"pdat\">p</span>/<span class=\"pdat\">pthresh</span>), </div> <p style=\"margin-top:0\">  where <span class=\"pdat\">pthresh</span> is 0.0005 by default.  The total score of a match is the sum of the p-scores of the hits making up the match.  <br><br>MCAST searches for all possible matches between the query motifs and the sequences in the database, and reports the matches with the largest scores in decreasing order. Three types of statistical confidence estimates (<i>p</i>-value, <i>E</i>-value, and <i>q</i>-value) are estimated for each score, and the reported matches can be filtered by applying <i>p</i>-value or <i>q</i>-value thresholds that you can specify.</footer>",

  "glam2scan": "<dfn>GLAM2Scan</dfn> searches sequences for matches to gapped <em>DNA</em> or <em>protein</em> GLAM2 motifs<span class=\"sample_output\"> (<a href=\"../doc/examples/glam2scan_example_output_files/glam2scan.html\">sample output</a> for <a href='http://meme-suite.org/meme-software/example-datasets/At.glam2'>GLAM2 motif</a> and <a href='http://meme-suite.org/meme-software/example-datasets/At.faa'>sequences</a>)</span>.<span class='manual_link'> See this <a href=\"../doc/glam2scan.html?man_type=web\">Manual</a> for more information.</span>  <footer><br><br>GLAM2SCAN finds matches, in a sequence database, to motifs discovered by GLAM2. Each match receives a score, indicating how well it fits the motif.  The GLAM2 motifs may be in either plain text (<code>glam2.txt</code>) or HTML (<code>glam2.html</code>) format.</footer>",

  "tomtom": "<dfn>Tomtom</dfn> compares one or more motifs against a database of known motifs (e.g., JASPAR). <dfn>Tomtom</dfn> will rank the motifs in the database and produce an alignment for each significant match<span class=\"sample_output\"> (<a href=\"../doc/examples/tomtom_example_output_files/tomtom.html\">sample output</a> for <a href='http://meme-suite.org/meme-software/example-datasets/STRGGTCAN.meme'>motif</a> and JASPAR CORE 2014 database)</span>.<span class='manual_link'> See this <a href=\"../doc/tomtom.html?man_type=web\">Manual</a> for more information.</span>  <footer><br><br>Tomtom searches one or more query motifs against one or more databases of target motifs (and their reverse complements when applicable), and reports for each query a list of target motifs, ranked by <i>p</i>-value.  The <i>E</i>-value and the <i>q</i>-value of each match is also reported.  The <i>q</i>-value is the minimal false discovery rate at which the observed similarity would be deemed significant.  The output contains results for each query, in the order that the queries appear in the input file.  <br><br>The web version of Tomtom provides a large number of target motif databases, or you can upload your own motifs in <a href=meme-format.html>MEME Motif Format</a>, or enter them by typing in many additional formats.  <br><br>For a given pair of motifs, the program considers all offsets between the motifs, while requiring a minimum number of overlapping positions. For a given offset, each overlapping position is scored using one of seven column similarity functions defined below. Columns in the query motif that don't overlap the target motif are assigned a score equal to the median score of the set of random matches to that column. <br><br>In order to compute the scores, Tomtom needs to know the frequencies of the letters of the sequence alphabet in the database being searched (the 'background' letter frequencies). By default, the background letter frequencies included in the query motif file are used. The scores of columns that overlap for a given offset are summed. This summed score is then converted to a <i>p</i>-value. The reported <i>p</i>-value is the minimal <i>p</i>-value over all possible offsets.  To compensate for multiple testing, each reported <i>p</i>-value is converted to an <i>E</i>-value by multiplying it by twice the number of target motifs.  As a second type of multiple-testing correction, <i>q</i>-values for each match are computed from the set of <i>p</i>-values and reported.</footer>",
  "annotated_sequences": "The MEME Suite <dfn>Motif Scanning</dfn> tools output their results as interactive HTML files. </br></br>You can see an example output from each of the motif scanning tools under <dfn>Motif Scanning</dfn> in the <dfn>Sample Outputs</dfn> menu on the left.",
  "aligned_motifs": "The MEME Suite <dfn>Tomtom</dfn> tool outputs its results as an interactive HTML file. </br></br>You can see an example <dfn>Tomtom</dfn> output under <dfn>Motif Comparison</dfn> in the <dfn>Sample Outputs</dfn> menu on the left, or by <b>clicking</b> here."
};

/*
 * Three types of URLS
 * 1) Relative URLs -
 *      Displayed all the time. These will have the path prefix appended. This
 *      type of URL is the default.
 * 2) Webserver URLs - 
 *      Only displayed when running on a webserver or if a webserver URL has
 *      been provided. When running on a webserver these should behave like
 *      relative URLs and have the path prefix appended, otherwise they should
 *      have the full server_url appended.
 * 3) Absolute URLs -
 *      Displayed all the time. These always go to the same site.
 *
 */
MemeMenu.data = {
  "header": {
    "title": "MEME Suite @VERSION@",
    "doc_url": "doc/overview.html",
    "web_url": "index.html"
  },
  "topics": [
    {
      "title": "Motif Discovery",
      "info": MemeMenu.blurbs["motif_discovery"],
      "server": true,
      "topics": [
        {
          "title": "MEME",
          "info": MemeMenu.blurbs["meme"],
          "url": "tools/meme",
          "server": true
        }, {
          "title": "DREME",
          "info": MemeMenu.blurbs["dreme"],
          "url": "tools/dreme",
          "server": true
        }, {
          "title": "MEME-ChIP",
          "info": MemeMenu.blurbs["memechip"],
          "url": "tools/meme-chip",
          "server": true
        }, {
          "title": "GLAM2",
          "info": MemeMenu.blurbs["glam2"],
          "url": "tools/glam2",
          "server": true
        }, {
          "title": "MoMo",
          "info": MemeMenu.blurbs["momo"],
          "url": "tools/momo",
          "server": true
        }
      ]
    }, {
      "title": "Motif Enrichment",
      "info": MemeMenu.blurbs["motif_enrichment"],
      "server": true,
      "topics": [
        {
          "title": "CentriMo",
          "info": MemeMenu.blurbs["centrimo"],
          "url": "tools/centrimo",
          "server": true
        }, {
          "title": "AME",
          "info": MemeMenu.blurbs["ame"],
          "url": "tools/ame",
          "server": true
        }, {
          "title": "SpaMo",
          "info": MemeMenu.blurbs["spamo"],
          "url": "tools/spamo",
          "server": true
        }, {
          "title": "GOMo",
          "info": MemeMenu.blurbs["gomo"],
          "url": "tools/gomo",
          "server": true
        }
      ]
    }, {
      "title": "Motif Scanning",
      "info": MemeMenu.blurbs["motif_scanning"],
      "server": true,
      "topics": [
        {
          "title": "FIMO",
          "info": MemeMenu.blurbs["fimo"],
          "url": "tools/fimo",
          "server": true
        }, {
          "title": "MAST",
          "info": MemeMenu.blurbs["mast"],
          "url": "tools/mast",
          "server": true
        }, {
          "title": "MCAST",
          "info": MemeMenu.blurbs["mcast"],
          "url": "tools/mcast",
          "server": true
        }, {
          "title": "GLAM2Scan",
          "info": MemeMenu.blurbs["glam2scan"],
          "url": "tools/glam2scan",
          "server": true
        }
      ]
    }, {
      "title": "Motif Comparison",
      "info": MemeMenu.blurbs["motif_comparison"],
      "server": true,
      "topics": [
        {
          "title": "Tomtom",
          "info": MemeMenu.blurbs["tomtom"],
          "url": "tools/tomtom",
          "server": true
        }
      ]
    }, {
      "title": "Manual",
      "info": "Read help pages for all the tools in the MEME Suite. This documentation can help you understand how the tools work and how to best use them.  Each help page describes the tool when run via the web, and also describes the inputs and outputs of the tool when you run it from the command line after installing MEME Suite on your computer.",
      "topics": [ 
        {
          "title": "OVERVIEW",
          "info": "See help pages for all the tools in the MEME Suite.",
          "url": "doc/overview.html"
        }, {
          "divider": true,
          "title": "Motif Discovery",
        }, {
          "title": "MEME",
          "info": "Read documentation on using MEME.",
          "url": "doc/meme.html",
        }, {
          "title": "DREME",
          "info": "Read documentation on using DREME.",
          "url": "doc/dreme.html"
        }, {
          "title": "MEME-ChIP",
          "info": "Read documentation on using MEME-ChIP.",
          "url": "doc/meme-chip.html"
        }, {
          "title": "GLAM2",
          "info": "Read documentation on using GLAM2.",
          "url": "doc/glam2.html"
        }, {
          "title": "MoMo",
          "info": "Read documentation on using MoMo.",
          "url": "doc/momo.html"
        }, {
          "divider": true,
          "title": "Motif Enrichment"
        }, {
          "title": "CentriMo",
          "info": "Read documentation on using CentriMo.",
          "url": "doc/centrimo.html"
        }, {
          "title": "AME",
          "info": "Read documentation on using AME.",
          "url": "doc/ame.html"
        }, {
          "title": "SpaMo",
          "info": "Read documentation on using SpaMo.",
          "url": "doc/spamo.html"
        }, {
          "title": "GOMo",
          "info": "Read documentation on using GOMo.",
          "url": "doc/gomo.html"
        }, {
          "divider": true,
          "title": "Motif Scanning"
        }, {
          "title": "FIMO",
          "info": "Read documentation on using FIMO.",
          "url": "doc/fimo.html"
        }, {
          "title": "MAST",
          "info": "Read documentation on using MAST.",
          "url": "doc/mast.html"
        }, {
          "title": "MCAST",
          "info": "Read documentation on using MCAST.",
          "url": "doc/mcast.html"
        }, {
          "title": "GLAM2Scan",
          "info": "Read documentation on using GLAM2Scan.",
          "url": "doc/glam2scan.html"
        }, {
          "divider": true,
          "title": "Motif Comparison"
        }, {
          "title": "Tomtom",
          "info": "Read documentation on using Tomtom.",
          "url": "doc/tomtom.html"
        }
      ]
    }, {
      "title": "Guides & Tutorials",
      "info": "Find out more about how to use the MEME Suite.",
      "topics": [
        {
          "title": "AME",
          "info": "Learn how to use AME to perform motif enrichment analysis on a set of sequences.",
          "url": "doc/ame-tutorial.html"
        }, {
          "title": "CentriMo",
          "info": "Learn how to use CentriMo to perform positional motif enrichment analysis on a set of sequences.",
          "url": "doc/centrimo-tutorial.html"
        }, {
          "title": "DREME",
          "info": "Learn how to use DREME to find short nucleotide motifs.",
          "url": "doc/dreme-tutorial.html"
          //TODO write article
        }, {
          "title": "FIMO",
          "info": "Learn how to use FIMO to find individual matches to each of your motifs in a set of sequences.",
          "url": "doc/fimo-tutorial.html"
        }, {
          "title": "GLAM2",
          "info": "Learn how to use GLAM2 to find GAPPED motifs.",
          "url": "doc/glam2_tut.html"
        }, {
          "title": "GT-Scan",
          "info": "Learn how to use GT-Scan to identify good targets for CRISPR/Cas and other genome targeting technologies.",
          "url": "http://bioinformatics.csiro.au/gt-scan/docs/manual",
          "absolute": true
        }, {
          "title": "MEME-ChIP",
          "info": "Learn how to analyze large DNA or RNA datasets using MEME-ChIP. It can analyze peak regions identified by ChIP-seq, cross-linking sites identified by CLIP-seq and related assays, as well as sets of genomic regions selected using other criteria. MEME-ChIP performs de novo motif discovery, motif enrichment analysis, motif location analysis and motif clustering, providing a comprehensive picture of the DNA or RNA motifs that are enriched in the input sequences.",
          "url": "http://www.nature.com/nprot/journal/v9/n6/full/nprot.2014.083.html",
          "absolute": true
        }/*, {
          "title": "Scan a genome",
          "info": "This article will tell you how to scan genome-sized sequence sets."
          //TODO write article
        }*/
      ]
    }, {
      "title": "Sample Outputs",
      "info": "See samples of the output produced by MEME Suite programs.",
      "topics": [
        {
          "divider": true,
          "title": "Motif Discovery",
        }, {
          "title": "MEME Sample",
          "info": "See sample output from MEME.",
          "url": "doc/examples/meme_example_output_files/meme.html"
        }, {
          "title": "DREME Sample",
          "info": "See sample output from DREME.",
          "url": "doc/examples/dreme_example_output_files/dreme.html"
        }, {
          "title": "MEME-ChIP Sample",
          "info": "See sample output from MEME-ChIP.",
          "url": "doc/examples/memechip_example_output_files/meme-chip.html"
        }, {
          "title": "GLAM2 Sample",
          "info": "See sample output from GLAM2.",
          "url": "doc/examples/glam2_example_output_files/glam2.html"
        }, {
          "title": "MoMo Sample",
          "info": "See sample output from MoMo.",
          "url": "doc/examples/momo_example_output_files/momo.html"
        }, {
          "divider": true,
          "title": "Motif Enrichment"
        }, {
          "title": "CentriMo Sample",
          "info": "See sample output from CentriMo.",
          "url": "doc/examples/centrimo_example_output_files/centrimo.html"
        }, {
          "title": "AME Sample",
          "info": "See sample output from AME.",
          "url": "doc/examples/ame_example_output_files/ame.html"
        }, {
          "title": "SpaMo Sample",
          "info": "See sample output from SpaMo",
          "url": "doc/examples/spamo_example_output_files/spamo.html"
        }, {
          "title": "GOMo Sample",
          "info": "See sample output from GOMo.",
          "url": "doc/examples/gomo_example_output_files/gomo.html"
        }, {
          "divider": true,
          "title": "Motif Scanning"
        }, {
          "title": "FIMO Sample",
          "info": "See sample output from FIMO.",
          "url": "doc/examples/fimo_example_output_files/fimo.html"
        }, {
          "title": "MAST Sample",
          "info": "See sample output from MAST.",
          "url": "doc/examples/mast_example_output_files/mast.html"
        }, {
          "title": "MCAST Sample",
          "info": "See sample output files from MCAST.",
          "url": "doc/examples/mcast_example_output_files/mcast.html"
        }, {
          "title": "GLAM2Scan Sample",
          "info": "See sample output files from GLAM2Scan.",
          "url": "doc/examples/glam2scan_example_output_files/glam2scan.html"
        }, {
          "divider": true,
          "title": "Motif Comparison"
        }, {
          "title": "Tomtom Sample",
          "info": "See sample output files from Tomtom.  ",
          "url": "doc/examples/tomtom_example_output_files/tomtom.html"
        }
      ]
    }, {
      "title": "File Format Reference",
      "info": "Read about the file formats that the MEME Suite uses.",
      "topics": [
        {
          "title": "FASTA Sequence",
          "info": "FASTA sequence format is the main sequence format read by tools in the MEME Suite.",
          "url": "doc/fasta-format.html"
        }, {
          "title": "MEME Motif",
          "info": "MEME motif format is the main motif format read by tools in the MEME Suite.",
          "url": "doc/meme-format.html"
        }, {
          "title": "Peptide-Spectra Match",
          "info": "Peptides identified from tandem mass spectra.",
          "url": "doc/psm-format.html"
        }, {
          "title": "Custom Alphabet",
          "info": "A custom alphabet can be specified to MEME Suite tools using this custom alphabet definition format.",
          "url": "doc/alphabet-format.html"
        }, {
          "title": "Markov Background Model",
          "info": "Background letter frequencies can be specified to MEME Suite tools using this Markov model based format.",
          "url": "doc/bfile-format.html"
        }, {
          "title": "Position-specific prior (PSP)",
          "info": "Priors (weights) on each position in each input sequence can be specified to MEME by a file in PSP format. These priors allow the user to bias the search for motifs by MEME.",
          "url": "doc/psp-format.html"
        }, {
          "title": "Dirichlet Mixtures",
          "info": "Dirichlet mixture priors capture the tendency of certain letters (e.g., protein residues) to align with each other. These priors allow the user to bias the search for motifs by MEME and GLAM2.",
          "url": "doc/dmix-format.html"
        }, {
          "title": "Other Supported Formats",
          "info": "See the complete list of formats that are documented in the MEME Suite.",
          "url": "doc/overview.html#formats"
        }
      ]
    }, {
      "title": "Databases",
      "info": "View information about the databases installed on this server.",
      "server": true,
      "topics": [
          {
            "title": "Motif Databases",
            "info": "View information on the motif databases installed on this server.",
            "server": true,
            "url": "db/motifs"
          }, {
            "title": "Sequence Databases",
            "info": "View information on the sequence databases installed on this server.",
            "server": true,
            "url": "db/sequences"
          }, {
            "title": "GOMo Sequence Databases",
            "info": "View information on the GOMo sequence databases installed on this server.",
            "server": true,
            "url": "db/gomo"
          }
        ]
    }, {
      "title": "Download & Install",
      "info": "Download a copy of the MEME Suite software and databases for running on your own computer.",
      "topics": [
        {
          "title": "Download MEME Suite and Databases",
          "info": "Download a copy of the MEME Suite for local installation, as well as motif databases and GOMo databases which can be used with it.",
          "url": "doc/download.html"
        }, {
          "title": "Copyright Notice",
          "info": "Read the copyright notice.",
          "url": "doc/copyright.html"
        }, {
          "title": "Installation Guide",
          "info": "Learn how to install the MEME Suite.",
          "url": "doc/install.html"
        }, {
          "title": "Release Notes",
          "info": "Find out what was added to the MEME Suite in this and previous releases.",
          "url": "doc/release-notes.html"
        }, {
          "title": "Release Announcement Group",
          "info": "Join this GOOGLE group to be notified of new MEME Suite releases and patches.",
          "url": "https://groups.google.com/forum/#!forum/meme-suite-releases",
          "absolute": true
        }, {
          "title": "Commercial Licensing",
          "info": "Acquire a licence for commercial use.",
          "url": "http://techtransfer.universityofcalifornia.edu/NCD/20911.html",
          "absolute": true
        }
      ]
    }, {
      "title": "Help",
      "info": "Get help using the MEME Suite.",
      "topics": [
        /*
        {
          "title": "Frequently Asked Questions",
          "info": "A list of frequently asked questions and their answers.",
          "url": "doc/general-faq.html"
          // This has been replaced by the Q&A forum
        },
        */
        {
          "title": "Q&A Forum",
          "info": "Ask questions on the forums.",
          "url": "https://groups.google.com/forum/#!forum/meme-suite",
          "absolute": true
        }, {
          "title": "Email Webmaster",
          "info": "Send the webmaster an email. This link will open your email client",
          "url": "mailto:@CONTACT@",
          "show_test": function(is_server, subs) {return subs["CONTACT"] !== ""},
          "absolute": true
        }, {
          "title": "Email Developers",
          "info": "Send the developers an email. This link will open your email client",
          "url": "mailto:@DEV_CONTACT@",
          "absolute": true
        }
      ]
    }, {
      "title": "Alternate Servers",
      "info": "View the current load at and submit jobs to the available MEME Suite servers.",
      "server": true,
      "topics": [
        {
          "title": "Main Server",
          "info": "Go to the main MEME Suite site.",
          "url": "http://meme-suite.org",
          "absolute": true
        },
        {
          "title": "->View Current Load",
          "info": "View the current load at the main MEME Suite site.",
          "url": "http://meme-suite.org/opal2/dashboard?command=statistics",
          "absolute": true
        },
        {
          "title": "Alternate Server",
          "info": "Go to the alternate MEME Suite site.",
          "url": "http://alternate.meme-suite.org",
          "absolute": true
        },
        {
          "title": "->View Current Load",
          "info": "View the current load at the alternate MEME Suite site.",
          "url": "http://alternate.meme-suite.org/opal2/dashboard?command=statistics",
          "absolute": true
        },
        {
          "title": "GenQuest (France)",
          "url": "http://tools.genouest.org/tools/meme/",
          "absolute": true
        }
      ]
    }, {
      "title": "Authors & Citing",
      "info": "Find out who contributed to the MEME Suite and how to cite the individual tools in publications.",
      "topics": [
        {
          "title": "Authors",
          "info": "The authors of the MEME Suite.",
          "url": "doc/authors.html"
        }, {
          "title": "Citing the MEME Suite",
          "info": "How to cite the MEME Suite or individual tools from the MEME Suite.",
          "url": "doc/cite.html"
        }
      ]
    }
  ],
  "opal_site": "http://sourceforge.net/projects/opaltoolkit/"
};
MemeMenu.script_loaded();
