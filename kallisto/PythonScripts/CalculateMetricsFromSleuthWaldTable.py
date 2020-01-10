import argparse
from sets import Set
import glob
from decimal import *
from sklearn.metrics import roc_auc_score

def ParseGoldStandard(goldfile,fdr):
    gold_dict = {}
    classifier = {}
    gold_open = open(goldfile,'r')
    gold_fields = gold_open.readline().strip().split()
    for line in gold_open:
        linelist =  line.strip().split()
        line_dict = dict(zip(gold_fields,linelist))
        if line_dict['qval'] != 'NA':
            gold_dict[line_dict['target_id']] = float(line_dict['qval'])
            if float(line_dict['qval']) <= fdr:
                classifier[line_dict['target_id']] = 1
            else:
                classifier[line_dict['target_id']] = 0
        elif line_dict['qval']  == 'NA':
            gold_dict[line_dict['target_id']] = 999
        else:
            raise Exception('unspecified condition')
    
    return gold_dict,classifier

def CalcAUC(classifier_dict,response_dict):
    classifier = []
    response = []
    for target in classifier_dict:
        classifier.append(classifier_dict[target])
        if target in response_dict:
            response.append(response_dict[target])
        else:
            resonse.append(0)
    auc = roc_auc_score(classifier,response)
    return auc

def ParseWald(truth_dict,evalfile,fdr):
    evalopen = open(evalfile,'r')
    fields = evalopen.readline().strip().split()
    true_positives = 0 ; false_positives = 0
    true_negatives = 0 ; false_negatives = 0
    response_dict = {}

    for line in evalopen:
        linelist = line.strip().split()
        line_dict = dict(zip(fields,linelist))
        
        if truth_dict[line_dict['target_id']] == 999:
            pass 
        else:
            if line_dict['qval'] == 'NA':
                line_dict['qval'] = 1
                response_dict[line_dict['target_id']] = 0
            else:
                ### response scores for auc calculations ###
                response_dict[line_dict['target_id']] = 1 - float(line_dict['qval'])

            ### performance metric tallying ###
            if truth_dict[line_dict['target_id']] > fdr and float(line_dict['qval']) > fdr:
                true_negatives +=1
            elif truth_dict[line_dict['target_id']] <= fdr and float(line_dict['qval']) > fdr:
                false_negatives +=1
            elif truth_dict[line_dict['target_id']] <=fdr and float(line_dict['qval']) <= fdr:
                true_positives +=1
            elif truth_dict[line_dict['target_id']] > fdr and float(line_dict['qval']) <= fdr:
                false_positives +=1
            else:
                raise Exception('invalid combination of fdr values')

    metric_dict = {}
    metric_dict['fp'] = false_positives/float(false_positives + true_negatives)
    metric_dict['fn'] = false_negatives/float(false_negatives + true_positives)
    metric_dict['sensitivity'] = true_positives/float(true_positives + false_negatives)
    metric_dict['specificity'] = true_negatives/float(true_negatives + false_positives)
    metric_dict['precision'] = true_positives/float(true_positives + false_positives) # aka ppv 
    
    return metric_dict,response_dict


if __name__=="__main__":

    parser = argparse.ArgumentParser(description='calculate spearman rank correlations among tpms for trimming study')
    parser.add_argument('-fdr','--false_discovery_rate',dest='fdr',type=float,help='fdr threshold')
    parser.add_argument('-sra','--sra-accession',dest='sra',type=str,help='sra accession id')
    parser.add_argument('-decount','--true-de-count',dest='de',type=int,help='min number of de transcripts from 2x125')
    parser.add_argument('-gene','--gene-level',dest='gene',action='store_true',help='specifies gene level wald tests to evaluate')
    opts = parser.parse_args()

    if opts.gene == True:
        condition_pairs = Set('_'.join(i.split('_')[3:-1]) for i in glob.glob('gene*%s*wald*tsv' % opts.sra) if 'summary' not in i)
    else:
        condition_pairs = Set('_'.join(i.split('_')[2:-1]) for i in glob.glob('%s*wald*tsv' % opts.sra) if 'summary' not in i)
    print 'condition pairs', condition_pairs
    if opts.gene == True:
        fout = open('%s_kallisto_waldgene_summary_fdr%s.tsv' % (opts.sra,opts.fdr),'w')
    else:
        fout = open('%s_kallisto_waldtranscript_summary_fdr%s.tsv' % (opts.sra,opts.fdr),'w')
    fout.write('accession\tcondition_pair\tstrategy\tfp\tfn\tsensitivity\tspecificity\tprecision\tfdr\tauc\n')
    
    for condition_pair in condition_pairs:
        if opts.gene == True:
            walds = glob.glob('gene*%s_[1,2]*x*[40,75,125]_%s_sleuthwald.tsv' % (opts.sra,condition_pair))
        else:
            walds = glob.glob('%s_[1,2]*x*[40,75,125]_%s_sleuthwald.tsv' % (opts.sra,condition_pair))
        print 'walds',walds
        gold_standard = [i for i in walds if '2x125' in i]
        print 'condition pair',condition_pair
        print 'gold standard',gold_standard
        if len(gold_standard) > 1:
            raise Exception('unspecified condition')

        truede = 0
        truth_dict,classifier_dict = ParseGoldStandard(gold_standard[0],opts.fdr)
 
        for feature in truth_dict:
            if truth_dict[feature] <= opts.fdr:
                truede +=1
        
        print 'truede = ', truede
        if truede >= opts.de:
            pe40 = [i for i in walds if '2x40' in i][0]
            se125 = [i for i in walds if '1x125' in i][0]
            se75 = [i for i in walds if '1x75' in i][0]

            pe40_metrics_dict,pe40_response_dict = ParseWald(truth_dict,pe40,opts.fdr)
            pe40auc = CalcAUC(classifier_dict,pe40_response_dict)
            fout.write('%s\t%s\tpe40\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (opts.sra,condition_pair,pe40_metrics_dict['fp'],pe40_metrics_dict['fn'],pe40_metrics_dict['sensitivity'],pe40_metrics_dict['specificity'],pe40_metrics_dict['precision'],opts.fdr,pe40auc)
)
            se75_metrics_dict,se75_response_dict = ParseWald(truth_dict,se75,opts.fdr)
            se75auc = CalcAUC(classifier_dict,se75_response_dict)
            fout.write('%s\t%s\tse75\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (opts.sra,condition_pair,se75_metrics_dict['fp'],se75_metrics_dict['fn'],se75_metrics_dict['sensitivity'],se75_metrics_dict['specificity'],se75_metrics_dict['precision'],opts.fdr,se75auc))

            se125_metrics_dict,se125_response_dict = ParseWald(truth_dict,se125,opts.fdr)
            se125auc = CalcAUC(classifier_dict,se125_response_dict) 
            fout.write('%s\t%s\tse125\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (opts.sra,condition_pair,se125_metrics_dict['fp'],se125_metrics_dict['fn'],se125_metrics_dict['sensitivity'],se125_metrics_dict['specificity'],se125_metrics_dict['precision'],opts.fdr,se125auc))
    
    fout.close()
