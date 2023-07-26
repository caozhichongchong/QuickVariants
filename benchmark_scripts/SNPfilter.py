import glob
import os
from datetime import datetime
import argparse
import pandas as pd
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='.',
                      metavar='.')
# optional input genome
optional.add_argument("-cluster",
                      help="a reference genome to run, default is all genomes",
                      type=str, default='',
                      metavar='genome1')
################################################## Definition ########################################################
args = parser.parse_args()
################################################### Set up ########################################################
# Set up cutoff
min_qual_for_call = 40 #Remove sample*candidate that has lower than this quality
min_maf_for_call = .9 #Remove sample*candidate, 0.9 for bacterial 0.7 for covid, 0.8 for correcting bacterial
contig_end = 10 # avoiding SNPs at the end of contigs
if 'covid' in args.i:
    min_maf_for_call = 0.7
min_cov = 3 # at least 3 reads mapped to POS
if 'human' in args.i:
    min_cov = 10
bowtievcfsuffix = 'flt.snp.vcf'#'.flt.snp.vcf'
mappervcfsuffix = 'mapper1.*vcf'
bowtiemappersamtovcfsuffix = '.mappersamtovcf.vcf'
Indel_depth_included = False # including indel depth when compute maf
Indel_POS_excluded = False # excluding SNPs covering indel POS
compare_to_baseline = False # compare to 0 SNP cases
Middle_only = False # only using middle reads
Depth_check = False # output a depth summary
if 'noindel' in args.i:
    compare_to_baseline = True
################################################### Function ########################################################
def filter_snp(depthall,depthsnp, POS, contig_length):
    MAF = depthsnp/depthall
    if (contig_length == 0 or (POS > contig_end and POS < contig_length - contig_end))and (MAF >= min_maf_for_call and depthsnp >= min_cov):#* 2):
        return True
    return False

def load_baselinesnp(bowtievcf):
    baseline_chrpos = set()
    for lines in open(bowtievcf + '.final.txt', 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        baseline_chrpos.add('%s\t%s'%(lines_set[0],lines_set[1]))
    return baseline_chrpos

def load_bowtie(bowtievcf):
    try:
        f1 = open(bowtievcf + '.final.txt', 'r')
    except IOError:
        baseline_chrpos = set()
        if compare_to_baseline:
            # load baseline
            numsnp = os.path.split(bowtievcf)[-1].split('.SNP.fasta')[0].split('.')[-1]
            tool = os.path.split(bowtievcf)[-1].split('.SNP.fasta.')[-1].split('.flt.snp.vcf')[0]
            if numsnp != '0':
                baselinegenome = bowtievcf.replace('%s.SNP.fasta'%(numsnp),'0.SNP.fasta')
                print('%s start loading baseline %s' % (datetime.now(), baselinegenome))
                if baselinegenome not in baseline_set:
                    baseline_set.setdefault(baselinegenome + tool,load_baselinesnp(baselinegenome))
                baseline_chrpos = baseline_set[baselinegenome + tool]
                print('%s finished loading baseline %s' % (datetime.now(), baselinegenome))
        # load depth of qualified indel to support snps on POSs of that indel
        indel_file = bowtievcf.replace('.flt.snp.vcf','.flt.indel.vcf')
        indel_depth = dict()
        if Indel_depth_included:
            print('%s start loading indel depth %s' % (datetime.now(), indel_file))
            for lines in open(indel_file, 'r'):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR, POS, USELESS, REF, ALT, QUAL = lines_set[:6]
                    if float(QUAL) >= min_qual_for_call:
                        POS = int(POS)
                        DPset = [float(x) for x in lines_set[9].split(':')[-1].replace('\n', '').split(',')]
                        # excluding DP of REF itself when considering indel depth
                        DPset[0] = 0
                        DP = sum(DPset)
                        for POSsub in range(POS, POS + len(REF)):
                            indel_depth.setdefault('%s\t%s'%(CHR,POSsub),DP)
        print('%s finished loading indel depth %s' % (datetime.now(), indel_file))
        # grep all snps
        vcfoutput = []
        snpoutput = ['CHR\tPOS\tREF\tALT\tDPall\tDPsnp\n']
        # load each line of snps
        snpline = 0
        for lines in open(bowtievcf, 'r'):
            snpline += 1
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, USELESS, REF, ALT, QUAL = lines_set[:6]
            ALT = ALT.split(',')
            if REF != 'N' and float(QUAL) >= min_qual_for_call:
                if baseline_chrpos== set() or '%s\t%s'%(CHR,POS) not in baseline_chrpos:
                    withsnp = False
                    DPset = [float(x) for x in lines_set[9].split(':')[-1].replace('\n', '').split(',')]
                    DP = sum(DPset) + indel_depth.get('%s\t%s'%(CHR,POS),0)
                    try:
                        contig_length = int(CHR.split('size')[1])
                    except IndexError:
                        try:
                            contig_length = int(CHR.split('_length_')[1].split('_')[0])
                        except IndexError:
                            contig_length = 0
                    for i in range(0, len(ALT)):
                        # each potential ALT
                        depthsnp = DPset[i+1]# skip REF
                        if filter_snp(DP, depthsnp, int(POS), contig_length):
                            # a qualified snp
                            withsnp = True
                            snpoutput.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (CHR, POS, REF, ALT[i], DP, depthsnp))
                    if withsnp:
                        vcfoutput.append(lines)
            if snpline % 1000000 == 0:
                print(datetime.now(), 'processed %s lines' % (snpline))
        # output qualified snps
        f1 = open(bowtievcf + '.final.vcf','w')
        f1.write(''.join(vcfoutput))
        f1.close()
        f1 = open(bowtievcf + '.final.txt', 'w')
        f1.write(''.join(snpoutput))
        f1.close()

def load_mapper(mappervcf):
    try:
        f1 = open(mappervcf + '.final.txt', 'r')
    except IOError:
        # grep all snps
        os.system('grep \';\' %s > %s.snp' % (mappervcf, mappervcf))
        vcfoutput = ['']
        snpoutput = ['CHR\tPOS\tREF\tALT\tDPall\tDPsnp\n']
        # load each line of snps
        snpline = 0
        baseline_chrpos = set()
        if compare_to_baseline:
            # load baseline
            numsnp = os.path.split(mappervcf)[-1].split('.SNP.fasta')[0].split('.')[-1]
            if numsnp != '0':
                baselinegenome = mappervcf.replace('%s.SNP.fasta' % (numsnp), '0.SNP.fasta')
                tool = os.path.split(mappervcf)[-1].split('.SNP.fasta.')[-1].split('.flt.snp.vcf')[0]
                print('%s start loading baseline %s' % (datetime.now(), baselinegenome))
                if baselinegenome not in baseline_set:
                    baseline_set.setdefault(baselinegenome + tool, load_baselinesnp(baselinegenome))
                baseline_chrpos = baseline_set[baselinegenome + tool]
                print('%s finished loading baseline %s' % (datetime.now(), baselinegenome))
        total_refSNP = os.path.basename(mappervcf).split('.SNP')[0].split('.')[-1]
        SNP = pd.read_csv(mappervcf + '.snp',sep='\t',header=None)
        SNP.columns = ['CHR','POS','REF','ALT','DP','DETAILS-MIDDLE','DETAILS-ENDS','SUPPORT']
        SNP = SNP[~SNP['REF'].isin(['N'])]
        SNPtrue = []
        totallines = SNP.shape[0]
        for indexline in range(0,totallines):
                snpline += 1
                CHR,POS,REF,ALT,DP,middleDP,endDP = list(SNP.iloc[indexline,:7])
                if POS > 0 and (not Indel_POS_excluded or SNP.iloc[min(indexline + 1, totallines - 1),1] > 0):
                    # not an indel, next POS not a start of an indel
                    ALT = ALT.split(',')
                    if baseline_chrpos == set() or '%s\t%s' % (CHR, POS) not in baseline_chrpos:
                            withsnp = False
                            DP = float(DP)
                            middleDP = [x.split(',') for x in middleDP.split(';')]
                            endDP = [x.split(',') for x in endDP.split(';')]
                            if Middle_only:
                                DP = sum([float(x[0]) + float(x[1]) for x in middleDP])
                            if DP > 0:
                                try:
                                    contig_length = int(CHR.split('size')[1])
                                except IndexError:
                                    try:
                                        contig_length = int(CHR.split('_length_')[1].split('_')[0])
                                    except IndexError:
                                        contig_length = 0
                                for i in range(0,len(ALT)):
                                    # each potential ALT
                                    if ALT[i] != '-':
                                        # not indel
                                        # using both middle and end depth
                                        if Middle_only:
                                            depthsnp = sum([float(x) for x in middleDP[i + 1]]) # skip REF
                                        else:
                                            depthsnp = sum([float(x) for x in middleDP[i + 1]]) + sum([float(x) for x in endDP[i + 1]]) # skip REF
                                        if Depth_check and depthsnp >= min_cov: #and depthsnp >= 0.5*DP:
                                            alloutput.append('%s\t%s\t%s\n'%(total_refSNP,depthsnp,DP))
                                        if filter_snp(DP, depthsnp, int(POS), contig_length):
                                            # a qualified snp
                                            withsnp = True
                                            snpoutput.append('%s\t%s\t%s\t%s\t%s\t%s\n'%(CHR,POS,REF,ALT[i],DP,depthsnp))
                                if withsnp:
                                    SNPtrue.append(indexline)
                    if snpline % 1000000 == 0:
                        print(datetime.now(), 'processed %s lines' % (snpline))
        # output qualified snps
        SNP.iloc[SNPtrue,:].to_csv(mappervcf + '.final.vcf',sep='\t',index = None)
        f1 = open(mappervcf + '.final.txt', 'w')
        f1.write(''.join(snpoutput))
        f1.close()

# find all snp vcfs
allvcf_bowtie = glob.glob('%s/SNP_model*/merge/%s*%s'%(args.i,args.cluster,bowtievcfsuffix))
allvcf_bowtie.sort()
allvcf_mapper = glob.glob('%s/SNP_model*/merge/%s*%s'%(args.i,args.cluster,mappervcfsuffix)) + \
                glob.glob('%s/SNP_model*/merge/%s*%s'%(args.i,args.cluster,bowtiemappersamtovcfsuffix))
allvcf_mapper = list(set(allvcf_mapper))
allvcf_mapper.sort()
# process bowtie vcfs
baseline_set = dict()
if Depth_check:
    alloutput = ['SNP\tDPsnp\tDP\n']
    fall = open('%s/SNP_model/alldepth.txt'%(args.i),'w')
for bowtievcf in allvcf_bowtie:
    if '.indel' not in bowtievcf:
        print('%s start processing %s' % (datetime.now(), bowtievcf))
        load_bowtie(bowtievcf)
        print('%s finished processing %s' % (datetime.now(), bowtievcf))
# process mapper vcfs
for mappervcf in allvcf_mapper:
    if '.indel' not in mappervcf and 'final' not in mappervcf:
        print('%s start processing %s' % (datetime.now(), mappervcf))
        if '.bcf.' in mappervcf:
            load_bowtie(mappervcf)
        else:
            load_mapper(mappervcf)
        print('%s finished processing %s' % (datetime.now(), mappervcf))
        if Depth_check:
            fall.write(''.join(alloutput))
            alloutput = []
if Depth_check:
    fall.close()
