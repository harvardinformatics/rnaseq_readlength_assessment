# RSEM
To estimate transcript and gene-level expression with RSEM, we first need to build an RSEM index. In doing so, we also specify which aligner we want to have RSEM use to generate alignments. For this study, we used (bowtie2)[http://bowtie-bio.sourceforge.net/bowtie2/index.shtml]. An example cmd for building the index, for <em>D. melanogaster</em> is as follows:

    rsem-prepare-reference -p 8 --bowtie2 --gtf $(pwd)/Drosophila_melanogaster.BDGP6.85.gtf $(pwd)/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa $(pwd)/Drosophila_melanogaster.BDGP6.dna_sm.toplevel

where -p indcates the number of threads. Note that both rsem and bowtie2 must be in your PATH for this command to work.

To estimate expression with paired-end reads, a generic example of our command line is as follows:

rsem-calculate-expression --bowtie2 -p 16 --paired-end R1.fastq R2.fastq /PATH/TO/Drosophila.melanogaster.BDGP6.dna.sm.toplevel <samplename>

where "<samplename>" indicates an outfile name, that is a concatenation of the sample name, and the sequencing strategy (so we don't overwrite results from other strategies!). As with kallisto, for single-end data, we need to specify library fragment mean and standard deviation, which we obtain from our kallisto quant analyses of 2x125. An example command line for single-end analysis with RSEM is as follows:

    rsem-calculate-expression --bowtie2 -p 16 --fragment-length-mean 300 --fragment-length-sd 60 R1.fastq /PATH/TO/Drosophila.melanogaster.BDGP6.dna.sm.toplevel <samplename>

 
