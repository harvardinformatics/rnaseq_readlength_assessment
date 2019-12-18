import argparse
from math import sqrt
from sets import Set
import glob
from sklearn.metrics import mean_squared_error

def ParseTsv(infile):
    fopen=open(infile,'r')
    tpm_dict= {}
    fopen.readline()
    for line in fopen:
        linelist=line.strip().split()
        tpm_dict[linelist[0]] = float(linelist[-1])
    
    return tpm_dict    



if __name__=="__main__":

    parser = argparse.ArgumentParser(description='calculate spearman rank correlations among tpms for trimming study')
    parser.add_argument('-sra','--accessionid',dest='sra',type=str,help='sra accession name')
    opts = parser.parse_args()
    
    strategies = ['2x40','1x125','1x75']
    quants = glob.glob('*2x125*quant')
    
    # get condition-by-replicate list
    samples = Set()
    for quant in quants:
        samples.add('_'.join(quant.split('_')[:-3]))
    samples = list(samples) 
    fout=open('kallisto_RMSE_by_strategy_%s.tsv' % opts.sra,'w')
    fout.write('accession\tsample\tstrategy\trmse\n')    

    for i in range(len(samples)):
        sample = samples[i]
        goldstdin = '%s_2x125_trimmed_quant/abundance.tsv' % sample
        goldtpm = ParseTsv(goldstdin)
        x = []
        for tsid in goldtpm:
            x.append(goldtpm[tsid])
            
        for i in range(3):
            tpmsample = '%s_%s_trimmed_quant/abundance.tsv' % (sample,strategies[i])
            tpmdict = ParseTsv(tpmsample)
            y = []
            for tsid in goldtpm:
                y.append(tpmdict[tsid])
            rmse = sqrt(mean_squared_error(x,y))
            
            fout.write('%s\t%s\t%s\t%s\n' % (opts.sra,sample,strategies[i],rmse))
     
 
    fout.close()

