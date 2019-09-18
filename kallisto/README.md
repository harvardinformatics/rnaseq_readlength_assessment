# kallisto

We estimated transcript-level expression with kallisto v. 0.43.1. First, for each reference genome used in our study, we extracted the reference transcript set using the RSEM<em> extract-reference-transcripts</em> utility. Next, we built a kallisto index for those transcripts, with an example for <em> D. melanogaster</em> as follows:

    kallisto index -i Drosophila_melanogaster.BDGP6.dna_sm.toplevel.transcripts.idx Drosophila_melanogaster.BDGP6.dna_sm.toplevel.transcripts.fa 

Expression was then estimated for each sample. An example, genetic command line for paired-end analysis is as follows:

    kallisto quant -i /PATH/TO/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.transcripts.idx -b 100 -t 6 -o <SAMPLENAME>_quant /PATH/TO/R1.fastq /PATH/TO/R2.fastq

Where the -b flag specifies the number of bootstraps (which are required for downstream DE analysis with sleuth), and -t indicated the number of threads. For single-end inference, one must specify the mean and standard deviation of the library fragment size. In our study, we use the paired-end 2x125 reads from each sample as our "truth", and obtain library size distribution information from the kallisto output for our truth set. In doing so, we round up to the nearest 50bp, given that, normally with single-end analysis, the true library size distribution is unknown. An example single-end command line is as follows:

    kallisto quant -i /PATH/TO/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.transcripts.idx --single --fragment-length=300 --sd=60 -b 100 -t 6 -o <SAMPLENAME>_quant /PATH/TO/R1.fastq

 
