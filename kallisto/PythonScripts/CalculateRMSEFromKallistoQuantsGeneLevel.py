import argparse
from math import sqrt
from sets import Set
import glob
#from scipy.stats import spearmanr
from sklearn.metrics import mean_squared_error

def ParseGeneTsv(infile):
    fopen=open(infile,'r')
    tpm_dict= {}
    fopen.readline()
    for line in fopen:
        linelist=line.strip().split()
        tpm_dict[linelist[1]] = float(linelist[0])
    
    return tpm_dict    



if __name__=="__main__":

    parser = argparse.ArgumentParser(description='calculate spearman rank correlations among tpms for trimming study')
    parser.add_argument('-sra','--accessionid',dest='sra',type=str,help='sra accession name')
    opts = parser.parse_args()
    
    strategies = ['2x40','1x125','1x75']
    quants = glob.glob('%s_genetpm_*2x125*tsv' % opts.sra)
    
    # get condition-by-replicate list
    samples = Set()
    for quant in quants:
        print quant
        samples.add('_'.join(quant.split('_')[:-1][2:]))
    samples = list(samples) 
    print samples
    fout=open('kallisto_geneRMSE_by_strategy_%s.tsv' % opts.sra,'w')
    fout.write('accession\tsample\tstrategy\trmse\n')    

    for i in range(len(samples)):
        sample = samples[i]
        goldstdin = '%s_genetpm_%s_2x125.tsv' % (opts.sra,sample)
        goldtpm = ParseGeneTsv(goldstdin)
        x = []
        for tsid in goldtpm:
            x.append(goldtpm[tsid])
            
        for i in range(3):
            tpmsample = '%s_genetpm_%s_%s.tsv' % (opts.sra,sample,strategies[i])
            tpmdict = ParseGeneTsv(tpmsample)
            y = []
            for tsid in goldtpm:
                y.append(tpmdict[tsid])
            rmse = sqrt(mean_squared_error(x,y)) 
            
            fout.write('%s\t%s\t%s\t%s\n' % (opts.sra,sample,strategies[i],rmse))
     
 
    fout.close()

