import argparse


def ParseTpmCorrelations(tpmcor):
    fopen = open(tpmcor,'r')
    fopen.readline()
    cor_dict = {}
    for line in fopen:
        linelist = line.strip().split('\t')
        accession,sample,strategy = linelist[0],linelist[1],linelist[2]
        cor_dict[(accession,sample,strategy)] = line.strip().split('\t')
    return cor_dict 

def ParseBowtieAlignmentStats(statsin):
    stats_dict = {}
    fopen = open(statsin,'r')
    fields = fopen.readline().strip().split('\t')
    for line in fopen:
        linedict = dict(zip(fields,line.strip().split('\t')))
        linedict['sample'] = '_'.join(linedict['sample'].replace('se_','').replace('SE_','').split('_')[:-1]) 
        if linedict['paired'] == 'y':
           linedict['strategy'] = '2x%s' % linedict['read_length'] 
        else:
            linedict['strategy'] = '1x%s' % linedict['read_length']    
        stats_dict[(linedict['accession'],linedict['sample'],linedict['strategy'])] = [linedict['nreads'],linedict['unaligned'],linedict['uniquely_aligned'],linedict['multi_aligned'],linedict['align_rate']]

    return stats_dict

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='calculate spearman rank correlations among tpms for trimming study')
    parser.add_argument('-tpm','--tpm-cor-data',dest='tpmcor',type=str,help='tpm correlation file')
    parser.add_argument('-alignstats','--bowtie2-alignment-stats',dest='alignstats',type=str,help='alignment stats summary file')
    opts = parser.parse_args()

    tpm_cor_dict =  ParseTpmCorrelations(opts.tpmcor)
    align_stats_dict = ParseBowtieAlignmentStats(opts.alignstats)

    fout = open('walignstats_%s' % opts.tpmcor,'w')
    fout.write('accession\tsample\tstrategy\tspearmanr\tpval\tnreads\tunaligned\tuniquely_aligned\tmulti_aligned\talign_rate\n')
    for entry in tpm_cor_dict:
        if entry in align_stats_dict:
            fout.write('%s\t%s\n' % ('\t'.join(tpm_cor_dict[entry]),'\t'.join(align_stats_dict[entry])))
        else:
            raise Exception('%s Align stats data not found' % str(entry))
fout.close()
