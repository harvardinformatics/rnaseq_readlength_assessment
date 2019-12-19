import argparse
from sets import Set
import glob
from math import sqrt
from sklearn.metrics import mean_squared_error

def rsem_parser(infile,gene=True):
    fopen = open(infile,'r')
    expression_dict = {}
    fields = fopen.readline().strip().split('\t')
    for line in fopen:
        linelist = line.strip().split('\t')
        linedict =  dict(zip(fields,linelist))
        if gene == True:
            expression_dict[linedict['gene_id']] = linedict['TPM']
        else:
            expression_dict[linedict['transcript_id']] = linedict['TPM']    
    return expression_dict

def CalcRMSE(x,y):
    if x.keys() != y.keys():
        RaiseException('results have different isoforms/genes')
    else:
        xvals = []
        yvals = []
        for key in x:
            xvals.append(float(x[key]))
            yvals.append(float(y[key]))
        rmse = sqrt(mean_squared_error(xvals,yvals))
    return rmse
        
 
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='calculate spearman rank correlations among tpms for trimming study')
    parser.add_argument('-sra','--accessionid',dest='sra',type=str,help='sra accession name')
    parser.add_argument('-gene','--gene-level',dest='gene',action='store_true',help='gene-level flag')
    opts = parser.parse_args()

    
    if opts.gene == True:
        level = 'genes'
    else:
        level = 'isoforms'    
   
    fout=open('%s_rsem_RMSE_by_strategy_%s.tsv' % (level,opts.sra),'w')
    fout.write('accession\tsample\tstrategy\trmse\n')
 
    pe_results = glob.glob('RSEM*[40,125]*%s*.results' % level)
    print pe_results 
    se_results = glob.glob('se*[125,75]*%s*.results' % level)
    
    if level == 'isoform':
        samples = samples=Set(['_'.join(i.replace('.isoforms.results','').split('_')[2:-1]) for i in pe_results])
    else:
        samples = samples=Set(['_'.join(i.replace('.genes.results','').split('_')[2:-1]) for i in pe_results])
    
    for sample in samples:
        print sample
        paired125 = glob.glob('RSEM*%s_125*%s.results' % (sample,level))
        #print paired125
        p125_dict = rsem_parser(paired125[0],gene=opts.gene)

        paired40 = glob.glob('RSEM*%s_40*%s.results' % (sample,level))    
        p40_dict = rsem_parser(paired40[0],gene=opts.gene)

        se125 = glob.glob('se*%s_125*%s.results' % (sample,level)) 
        s125_dict = rsem_parser(se125[0],gene=opts.gene)

        se75 = glob.glob('se*%s_75*%s.results' % (sample,level))
        s75_dict = rsem_parser(se75[0],gene=opts.gene)    

        p40_rmse= CalcRMSE(p125_dict,p40_dict)
        fout.write('%s\t%s\t2x40\t%s\n' % (opts.sra,sample,p40_rmse))

        s75_rmse = CalcRMSE(p125_dict,s75_dict)
        fout.write('%s\t%s\t1x75\t%s\n' % (opts.sra,sample,s75_rmse))        

        s125_rmse = CalcRMSE(p125_dict,s125_dict)
        fout.write('%s\t%s\t1x125\t%s\n' % (opts.sra,sample,s125_rmse))

    fout.close()
