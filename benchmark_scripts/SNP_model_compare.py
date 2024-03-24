################################################## END ########################################################
################################################### SET PATH ########################################################
# compare bowtie call SNPs VS mapper call SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import argparse
from datetime import datetime
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='Gut_microbiome_benchmark/',
                      metavar='.')
# optional input genome
optional.add_argument("-cluster",
                      help="a reference genome to run, default is all genomes",
                      type=str, default='',
                      metavar='genome1')
################################################## Definition ########################################################
args = parser.parse_args()
# set up path
bowtie_checkFP = True
mapper_checkFP = True
# function
def load_vcf_ref(vcf_file):
    vcf_input = set()
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            REF = lines_set[3] # mutation introduced in reference
            CHRPOS = '%s\t%s\t%s\t'%(CHR,POS,REF)
            vcf_input.add(CHRPOS)
    return vcf_input

def load_vcf(vcf_file):
    vcf_input = set()
    for lines in open(vcf_file.replace('.final.vcf','.final.txt')):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            REF = lines_set[2] # mutation introduced in reference
            CHRPOS = '%s\t%s\t%s\t' % (CHR, POS, REF)
            vcf_input.add(CHRPOS)
    return vcf_input

def compare_vcf(vcf,print_vcf_set=False):
    summary_file_output = []
    print(datetime.now(), 'started processing', vcf)
    vcf_input1 = load_vcf(vcf)
    if '.mapper1.vcf' in vcf:
        print_vcf_set = True
    FN_diff = compare_vcf_detail(vcf, vcf_inputref, vcf_input1,print_vcf_set, 'FN') # vcf diff in ref not in 1, FN
    FP_diff = compare_vcf_detail(vcf, vcf_input1, vcf_inputref,print_vcf_set, 'FP') # vcf diff in 1 not in ref, FP
    summary_file_output.append('%s\t%s\t%s\t%s\n'%(
        os.path.split(vcf)[-1],
        FP_diff,
        FN_diff,total_SNP_ref
    ))
    f1 = open(summary_file, 'a')
    f1.write(''.join(summary_file_output))
    f1.close()
    print(datetime.now(), 'finished processing', vcf)

def compare_vcf_detail(vcf, vcf_input1, vcf_input2, print_vcf_set, TAG):
    # CHRPOS in 1 not in 2
    vcf_set = ['\t'.join(CHRPOS.split('\t')[:2]) + '\t' for CHRPOS in vcf_input1 if CHRPOS not in vcf_input2]
    if bowtie_checkFP and 'mapper' in vcf:
        f1 = open(vcf + '.' + TAG, 'w')
        f1.write('\n'.join(vcf_set))
        f1.close()
    if print_vcf_set and '.bcf.' not in vcf:
        # print max 100 FPs/FNs
        temp_output = os.path.join(output_dir, 'grep.temp.txt')
        f1 = open(temp_output, 'w')
        f1.write('\n'.join(vcf_set))
        f1.close()
        os.system('grep -T -f %s %s --no-group-separator > %s' % (
            temp_output,
            vcf.replace('.final.vcf',''),
            vcf + '.temp'))
        os.system('sort -k3 -n %s | sort -k2 > %s' %
                  (vcf + '.temp', vcf + '.' + TAG)
                  )
        os.system('rm -rf %s' % (vcf + '.temp'))
        os.system('rm -rf %s' % (temp_output))
    return len(vcf_set)

allresultdir = glob.glob('%s/SNP_model*/'%(args.i))
for resultdir in allresultdir:
    resultdirfilename = os.path.split(resultdir)[-1]
    output_dir = '%s/merge/'%(resultdir)
    ref_dir = '%s/data/'%(resultdir)
    vcf_ref = glob.glob(ref_dir + '/%s*.snp.txt'%(args.cluster))
    # summarize results
    summary_file = output_dir + '/model.sum.%s.txt'%(resultdirfilename)
    f1=open(summary_file,'w')
    f1.write('sample\tFP\tFN\ttotal_refSNP\n')
    f1.close()
    vcf_ref.sort()
    print(vcf_ref)
    for vcf_ref_file in vcf_ref:
        vcf_ref_file_name = os.path.split(vcf_ref_file)[-1].split('.snp.txt')[0]
        print('process vcfs for reference %s'%(vcf_ref_file_name))
        # load ref
        vcf_inputref = load_vcf_ref(vcf_ref_file)
        total_SNP_ref = len(vcf_inputref)
        allvcf = glob.glob(output_dir + vcf_ref_file_name + '*.final.vcf')
        for vcf in allvcf:
            print('processing %s' % (vcf))
            compare_vcf(vcf)