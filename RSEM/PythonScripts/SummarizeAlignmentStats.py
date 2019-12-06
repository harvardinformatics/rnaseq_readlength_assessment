import sys
import glob

def ParseStderr(filename):
    fopen = open(filename,'r')
    nreads,pairing,unaligned,unique,multi,overall = fopen.read().split('\n')[:-1]
    nreads = nreads.split()[0] 
    unaligned = unaligned.split('(')[1].split(')')[0].replace('%','')
    unique = unique.split('(')[1].split(')')[0].replace('%','')
    multi = multi.split('(')[1].split(')')[0].replace('%','')
    overall = overall.split()[0].replace('%','')
    stats_dict = {'jobid' : filename.replace('.e','').split('_')[-1],'nreads' : nreads,'unaligned' : unaligned,'unique' : unique,'multi' : multi,'overall': overall}
    if '_se' in filename:
        stats_dict['paired-end'] = 'n'
    else:
        stats_dict['paired-end'] = 'y'
    return stats_dict
 

def ParseStdOut(filename):
    outopen = open(filename,'r')
    lines = outopen.readlines()[0:2]
    for line in lines:
        if 'bowtie2' in line:
            sample = line.strip().split()[-2].split('/')[1].replace('RSEM_bt2_','').replace('.bam','') 
    return sample

errlogs = glob.glob("*.e")
fout = open('%s_bowtie2_alignmentstats_summary.tsv' % sys.argv[1],'w')
fout.write('accession\tsample\tjobid\tread_length\tpaired\tnreads\tunaligned\tuniquely_aligned\tmulti_aligned\talign_rate\n')
for errlog in errlogs:
    err_stats =  ParseStderr(errlog)
    stdout = '%s.o' % errlog.split('.')[0]
    sample = ParseStdOut(stdout)
    fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (sys.argv[1],sample,err_stats['jobid'],sample.split('_')[-1],err_stats['paired-end'],err_stats['nreads'],err_stats['unaligned'],err_stats['unique'],err_stats['multi'],err_stats['overall']))

fout.close()
    
