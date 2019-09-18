# RSEM
To estimate transcript and gene-level expression with RSEM, we first need to build an RSEM index. In doing so, we also specify which aligner we want to have RSEM use to generate alignments. For this study, we used (bowtie2)[http://bowtie-bio.sourceforge.net/bowtie2/index.shtml]. An example cmd for building the index, for <em>D. melanogaster</em> is as follows:

    rsem-prepare-reference -p 8 --bowtie2 --gtf $(pwd)/Drosophila_melanogaster.BDGP6.85.gtf $(pwd)/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa $(pwd)/Drosophila_melanogaster.BDGP6.dna_sm.toplevel

where -p indcates the number of threads. Note that both rsem and bowtie2 must be in your PATH for this command to work.
