import argparse
from sets import Set
import glob
from scipy.stats import spearmanr

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

def CalcSpearman(x,y):
    if x.keys() != y.keys():
        RaiseException('results have different isoforms/genes')
    else:
        xvals = []
        yvals = []
        for key in x:
            xvals.append(float(x[key]))
            yvals.append(float(y[key]))
        cor,pval = spearmanr(xvals,yvals)
    return cor,pval
        
 
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='calculate spearman rank correlations among tpms for trimming study')
    parser.add_argument('-sra','--accessionid',dest='sra',type=str,help='sra accession name')
    parser.add_argument('-gene','--gene-level',dest='gene',action='store_true',help='gene-level flag')
    opts = parser.parse_args()

    
    if opts.gene == True:
        level = 'genes'
    else:
        level = 'isoforms'    
   
    fout=open('%s_rsem_tpmcorrelations_by_strategy_%s.tsv' % (level,opts.sra),'w')
    fout.write('%s\tsample\tstrategy\tspearmanr\tpval\n' % opts.sra)
 
    #strategies = ['2x40','1x125','1x75']
    
    pe_results = glob.glob('RSEM*[40,125]*%s*.results' % level)
    print 'pe results', pe_results
    se_results = glob.glob('se*[125,75]*%s*.results' % level)
    print 'se results',se_results
    if level == 'isoform':
        samples = samples=Set(['_'.join(i.replace('.isoforms.results','').split('_')[2:-1]) for i in pe_results])
    else:
        samples = samples=Set(['_'.join(i.replace('.genes.results','').split('_')[2:-1]) for i in pe_results])
    
    for sample in samples:
        print sample
        #KO_Inf_rep3_125
        #RSEM_bt2_KO_Inf_rep1_125_adapt.genes.results
        paired125 = glob.glob('RSEM*%s_125*%s.results' % (sample,level))
        print 'paired 125',paired125
        p125_dict = rsem_parser(paired125[0],gene=opts.gene)

        paired40 = glob.glob('RSEM*%s_40*%s.results' % (sample,level))    
        print 'paired 40',paired40
        p40_dict = rsem_parser(paired40[0],gene=opts.gene)

        se125 = glob.glob('se*%s_125*%s.results' % (sample,level)) 
        print 'se125',se125
        s125_dict = rsem_parser(se125[0],gene=opts.gene)

        se75 = glob.glob('se*%s_75*%s.results' % (sample,level))
        print 'se75', se75
        s75_dict = rsem_parser(se75[0],gene=opts.gene)    

        p40_cor,p40_p = CalcSpearman(p125_dict,p40_dict)
        fout.write('%s\t%s\t2x40\t%s\t%s\n' % (opts.sra,sample,p40_cor,p40_p))

        s75_cor,s75_p = CalcSpearman(p125_dict,s75_dict)
        fout.write('%s\t%s\t1x75\t%s\t%s\n' % (opts.sra,sample,s75_cor,s75_p))        

        s125_cor,s125_p = CalcSpearman(p125_dict,s125_dict)
        fout.write('%s\t%s\t1x125\t%s\t%s\n' % (opts.sra,sample,s125_cor,s125_p))

    fout.close()
