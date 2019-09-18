import argparse
from sets import Set
import glob
from decimal import *
#getcontext().prec=6


def ParseWald(truth_dict,evalfile,fdr):
    #observations = len(truth_dict)
    #positives = 0
    #negatives = 0
    #for target in truth_dict:
        #if float(truth_dict[target]) <= fdr:
            #positives +=1  
        #else:
            #negatives +=1

    evalopen = open(evalfile,'r')
    fields = evalopen.readline().strip().split()
    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
    for line in evalopen:
        linelist = line.strip().split()
        line_dict = dict(zip(fields,linelist))
        if line_dict['qval'] == 'NA' and truth_dict[line_dict['target_id']] == 999:
            pass # not expressed
        else:
            if line_dict['qval'] == 'NA':
                line_dict['qval'] = 1
    
            if truth_dict[line_dict['target_id']] > fdr and float(line_dict['qval']) > fdr:
                true_negatives +=1
            elif truth_dict[line_dict['target_id']] <= fdr and float(line_dict['qval']) > fdr:
                false_negatives +=1
            elif truth_dict[line_dict['target_id']] <=fdr and float(line_dict['qval']) <= fdr:
                true_positives +=1
            elif truth_dict[line_dict['target_id']] > fdr and float(line_dict['qval']) <= fdr:
                false_positives +=1
            else:
                print 'wtf'
    metric_dict = {}
    metric_dict['fp'] = false_positives/float(false_positives + true_negatives)
    metric_dict['fn'] = false_negatives/float(false_negatives + true_positives)
    metric_dict['sensitivity'] = true_positives/float(true_positives + false_negatives)
    metric_dict['specificity'] = true_negatives/float(true_negatives + false_positives)
    metric_dict['precision'] = true_positives/float(true_positives + false_positives) # aka ppv 
    return metric_dict


if __name__=="__main__":

    parser = argparse.ArgumentParser(description='calculate spearman rank correlations among tpms for trimming study')
    parser.add_argument('-fdr','--false_discovery_rate',dest='fdr',type=float,help='fdr threshold')
    parser.add_argument('-sra','--sra-accession',dest='sra',type=str,help='sra accession id')
    parser.add_argument('-decount','--true-de-count',dest='de',type=int,help='min number of de transcripts from 2x125')
    opts = parser.parse_args()

    condition_pairs = Set('_'.join(i.split('_')[3:-1]) for i in glob.glob('genelevel*%s*wald*tsv' % opts.sra) if 'summary' not in i)
    fout = open('%s_kallisto_waldgene_summary_fdr%s.tsv' % (opts.sra,opts.fdr),'w')
    fout.write('accession\tcondition_pair\tstrategy\tfp\tfn\tsensitivity\tspecificity\tprecision\tfdr\n')
    for condition_pair in condition_pairs:
        walds = glob.glob('genelevel_*%s*%s_sleuthwald*tsv' % (opts.sra,condition_pair))
        gold_standard = [i for i in walds if '2x125' in i] 
        truth_dict = {}
        gold_open = open(gold_standard[0],'r')
        gold_fields = gold_open.readline().strip().split()
        for line in gold_open:
            linelist =  line.strip().split()
            line_dict = dict(zip(gold_fields,linelist))
            if linelist[1] != 'NA':
                truth_dict[line_dict['target_id']] = float(line_dict['qval'])
            elif linelist[1] == 'NA':
                truth_dict[line_dict['target_id']] = 999
            else:
                RaiseException('unspecified condition')
        ### 
        truede = 0
        for feature in truth_dict:
            if truth_dict[feature] <= opts.fdr:
                truede +=1
        
        print 'truede = ', truede
        if truede >= opts.de:
        ###
            pe40 = [i for i in walds if '2x40' in i][0]
            se125 = [i for i in walds if '1x125' in i][0]
            se75 = [i for i in walds if '1x75' in i][0]

            pe40_metrics_dict = ParseWald(truth_dict,pe40,opts.fdr)
            fout.write('%s\t%s\tpe40\t%s\t%s\t%s\t%s\t%s\t%s\n' % (opts.sra,condition_pair,pe40_metrics_dict['fp'],pe40_metrics_dict['fn'],pe40_metrics_dict['sensitivity'],pe40_metrics_dict['specificity'],pe40_metrics_dict['precision'],opts.fdr))

            se75_metrics_dict = ParseWald(truth_dict,se75,opts.fdr)
            fout.write('%s\t%s\tse75\t%s\t%s\t%s\t%s\t%s\t%s\n' % (opts.sra,condition_pair,se75_metrics_dict['fp'],se75_metrics_dict['fn'],se75_metrics_dict['sensitivity'],se75_metrics_dict['specificity'],se75_metrics_dict['precision'],opts.fdr)
)
            se125_metrics_dict = ParseWald(truth_dict,se125,opts.fdr) 
            fout.write('%s\t%s\tse125\t%s\t%s\t%s\t%s\t%s\t%s\n' % (opts.sra,condition_pair,se125_metrics_dict['fp'],se125_metrics_dict['fn'],se125_metrics_dict['sensitivity'],se125_metrics_dict['specificity'],se125_metrics_dict['precision'],opts.fdr))
    fout.close()
