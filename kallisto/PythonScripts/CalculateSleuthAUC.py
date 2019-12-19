import argparse
from sets import Set
from sklearn.metrics import roc_auc_score

def ParseGoldStandard(goldresults,fdr):
    classifier_dict = {}
    fopen = open(goldresults,'r')
    true_fields = 'target_id num_aggregated_transcripts sum_mean_obs_counts pval qval coeff'.split()
    fields = fopen.readline().split()
    if fields != true_fields:
        raise Exception('column names do notmatch expected resuilts for sleuth output')    
    else:
        for line in fopen:
            linedict = dict(zip(fields,line.strip().split()))
            if linedict['qval'] != 'NA':
                if float(linedict['qval'])<=fdr:
                    classifier_dict[linedict['target_id']] = 1
                else:
                    classifier_dict[linedict['target_id']] = 0
    return classifier_dict        

def ParseTargetWaldTable(waldtable):
    score_dict = {}
    fopen = open(waldtable,'r')
    true_fields = 'target_id num_aggregated_transcripts sum_mean_obs_counts pval qval coeff'.split()
    fields = fopen.readline().split()
    if fields != true_fields:
        raise Exception('column names do notmatch expected resuilts for sleuth output')    
    else:
        for line in fopen:
            linedict = dict(zip(fields,line.strip().split()))
            if linedict['qval'] != 'NA':
                score_dict[linedict['target_id']] = 1 - float(linedict['qval'])  

    return score_dict  

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='calculate auc from sleuth wald test')
    parser.add_argument('-sra','--accessionid',dest='sra',type=str,help='sra accession name')
    parser.add_argument('-g','--2x125-wald-result-table',dest='waldgold',type=str,help='sleuth 2x125 output')
    parser.add_argument('-fdr','--gold-standard-fdr',dest='fdr',type=float,help='value for defining DE in target table')
    parser.add_argument('-w','--wald-target-table',dest='wald',type=str,help='sleuth wald test output')
    opts = parser.parse_args()
    filelist = opts.wald.split('_')
    test = filelist[-2]
    strategy = filelist[-3] 
    classifier_dict = ParseGoldStandard(opts.waldgold,opts.fdr)
    target_dict = ParseTargetWaldTable(opts.wald)
    
    classifier = []
    response = []
    for target in classifier_dict:
        classifier.append(classifier_dict[target])
        if target in target_dict:
            response.append(target_dict[target])
        else:
            response.append(0)

    auc = roc_auc_score(classifier,response)
    print opts.sra,test,strategy,auc

