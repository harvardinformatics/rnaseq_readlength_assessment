# RSEM
To estimate transcript and gene-level expression with RSEM, we first need to build an RSEM index. In doing so, we also specify which aligner we want to have RSEM use to generate alignments. For this study, we used (bowtie2)[http://bowtie-bio.sourceforge.net/bowtie2/index.shtml]. An example cmd for building the index, for <em>D. melanogaster</em> is as follows:

    rsem-prepare-reference -p 8 --bowtie2 --gtf $(pwd)/Drosophila_melanogaster.BDGP6.85.gtf $(pwd)/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa $(pwd)/Drosophila_melanogaster.BDGP6.dna_sm.toplevel

where -p indcates the number of threads. Note that both rsem and bowtie2 must be in your PATH for this command to work.

To estimate expression with paired-end reads, a generic example of our command line is as follows:

rsem-calculate-expression --bowtie2 -p 16 --paired-end R1.fastq R2.fastq /PATH/TO/Drosophila.melanogaster.BDGP6.dna.sm.toplevel <samplename>

where "<samplename>" indicates an outfile name, that is a concatenation of the sample name, and the sequencing strategy (so we don't overwrite results from other strategies!). As with kallisto, for single-end data, we need to specify library fragment mean and standard deviation, which we obtain from our kallisto quant analyses of 2x125. An example command line for single-end analysis with RSEM is as follows:

    rsem-calculate-expression --bowtie2 -p 16 --fragment-length-mean 300 --fragment-length-sd 60 R1.fastq /PATH/TO/Drosophila.melanogaster.BDGP6.dna.sm.toplevel <samplename>

## Differential expression with limma (voom)
For our truth set 2x125 paired end reads, as well our evaluated strategies (2x40,1x75,1x125), we carried out pairwise Wald tests of differential expression using limma voom. We created separate expression matrices for RSEM outputs at the transcript and gene levels using the RSEM <em>rsem-generate-data-matrix</em> tool. Then, depending upon the pairwise analysis we conducted, we subsetted this matrix to generate the appropriate design matrix. As with sleuth-based differential expression analyses, we do not provide here our python code for generating R scripts, as they depend upon the naming convention one uses in prior pipeline steps. The arguments passed to execute limma analyses are identical between gene and transcript levels, only varying in which expression matrix one uses. Thus we provide one generic example R script, [SRP096374_TSvsWT_2x40_limma.Rscript](https://github.com/harvardinformatics/rnaseq_readlength_assessment/blob/master/RSEM/Rscripts/SRP096374_TSvsWT_2x40_limma.Rscript). 

We then calcuate differential expression test performance metrics, using 2x125 as our truth set. We do this using [RSEM_CalculateMetricsLimmaOutput.py](https://github.com/harvardinformatics/rnaseq_readlength_assessment/blob/master/RSEM/PythonScripts/RSEM_CalculateMetricsLimmaOutput.py), which iterates over the Wald test tables generated for all the pairwise tests within an SRA accession. We can then concatenate these tables in order to facilitate visualizion of global, multi-accession patterns. Note, the arguments passed to the python <em>glob </em> function may need to be changed,depending upon the naming convention used for fastq files.
