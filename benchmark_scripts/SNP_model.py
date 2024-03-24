# step 1 modelling SNPs and test 2 methods
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of WGS data",
                      type=str, default='Gut_microbiome_benchmark/original_data/',
                      metavar='input/')
# optional output setup
optional.add_argument("-fa",
                      help="file extension of fasta files",
                      type=str, default='.fasta.corrected',
                      metavar='.fasta')
optional.add_argument("-fq",
                      help="file extension of fastq files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')

optional.add_argument("-s",
                      help="a folder for your mapper",
                      type=str, default='SNP_model_script/',
                      metavar='SNP_model_script/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='Gut_microbiome_benchmark/',
                      metavar='.')
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 5)",
                      metavar="1 or more", action='store', default=5, type=int)
optional.add_argument("-indel",
                      help="whether to insert indels",
                      type=str, default='True',
                      metavar='False or True')

################################################## Definition ########################################################
# setting match score as 0 for all
# mismatch 6, match 0 -> mismatch 4, match 2
penalty = [6,5,3,0,1] # bowtie #mismatch penalty, gap open penalty, extension penalty, match score (not used for end-to-end), N dna penalty
#penalty = [4+2,4,2,0,2-1] # minimap, match score = 2, np score -> change to N dna penalty, turning 2 into -2 = -4
#penalty = [4+1,6,1,0,0+1] # bwa, match score = 1, no N dna penalty
#penalty = [1,2,0.5,0,0.1] # mapper
args = parser.parse_args()
input_script = args.s
genome_root = args.i
output_dir = args.o + '/SNP_model'
genome_name = args.fa
fastq_name=args.fq
fastq_name2=args.fq.replace('1','2')
threadstouse=args.t
cause_SNP = False
mapping_file = True
time_evaluation = False
penalty_test = False
mutation_repeat = False
mappertovcfonly = False
correct_genome = False
run_bcf = True
parallel_test = False
if time_evaluation:
    input_script_sub = '%s/SNP_model_newtime'%(input_script)
elif parallel_test:
    input_script_sub = '%s/SNP_model_parallel'%(input_script)
else:
    input_script_sub = '%s/SNP_model_MGsimbig'%(input_script)
output_dir_temp = output_dir
if time_evaluation:
    output_dir_temp = '/dev/shm/temp/'
    try:
        os.mkdir(output_dir_temp)
    except IOError:
        pass
    try:
        os.mkdir(output_dir_temp+'/bwa')
    except IOError:
        pass
    try:
        os.mkdir(output_dir_temp+'/merge')
    except IOError:
        pass
latest_mapper = glob.glob('%s/mapper-1*.jar'%(args.s))
latest_mapper.sort()
latest_mapper=latest_mapper[-1]
latest_mappersamtovcf = glob.glob('%s/quick-variants*.jar'%(args.s))
latest_mappersamtovcf.sort()
latest_mappersamtovcfold = latest_mappersamtovcf[0]
latest_mappersamtovcf=latest_mappersamtovcf[-1]
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# Set up N or S
N_S_set = dict()
N_S_set['N']=0
N_S_set['S']=1
purines=['A','G']
pyrimidines=['C','T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
# Set up mutation rate

if correct_genome or 'covid' in output_dir:
    genome_name = genome_name.replace('.corrected','')
    mutation_repeat = True
try:
    os.mkdir(args.o)
except IOError:
    pass

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/data')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/merge')
except IOError:
    pass

os.system('rm -r %s'%(input_script_sub))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# function
def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq[position - 1]=ALT
    return seq

def translate(seq):
    seq = Seq(''.join(seq))
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        if ALT_frq > 0:
            ALT_set.setdefault(ALT_frq, set())
            ALT_set[ALT_frq].add(alleles)
            ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def transitions(REF,ALT):
    if REF in pyrimidines:
        REF = complement[REF]
        ALT = complement[ALT]
    return '%s-%s'%(REF,ALT)

def contig_to_gene(CHR, POS):
    return []

def loaddatabase(database_aa,database):
    # load database seq
    Length = []
    reference_database = os.path.split(database_aa)[-1]
    print('reference database_aa set as %s' % (reference_database))
    Ref_seq = dict()
    Input_seq = dict()
    Input_id = []
    corrected_length = dict()
    for lines in open(database + '.txt','r'):
        record_id = lines.split('\t')[0]
        corrected_length.setdefault(record_id,0)
        corrected_length[record_id] += 1
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        seq_length = len(record_seq)
        Input_seq.setdefault(record_id, record_seq)
        Input_id.append(record_id)
        Length.append(seq_length - corrected_length.get(record_id,0))
    for record in SeqIO.parse(database_aa, 'fasta'):
        record_id = str(record.id)
        record_seq = list(str(record.seq))
        Ref_seq.setdefault(record_id, record_seq)
    return [Ref_seq,Length,Input_seq,Input_id]

def modelindel(seq,Chr,indel_set):
    SNP_output = []
    indel_set.sort()
    for position in indel_set:
        total_length = len(seq)
        if position < total_length:
            REF = seq[position]
            gene_info = contig_to_gene(Chr, position)
            temp_ALT = ['A', 'T', 'G', 'C']
            indel_size = random.choices(indel_nonorf, k=1)[0]
            if indel_size > 0:# insertion on ref
                ALT = random.choices(temp_ALT, k=indel_size)
                seq = seq[:position + 1] + ALT + seq[position + 1:]
                # including the REF, add insertion after the REF
                temp_line = [Chr, str(position + 1), REF, REF + ''.join(ALT)]
            else:# deletion on ref
                REF = ''.join(seq[(position+indel_size):position])
                del seq[(position+indel_size):position]
                temp_line = [Chr, str(position + 1 +indel_size), REF, '-'*(-indel_size)]
            SNP_output.append('\t'.join(temp_line) + '\n')
        else:
            print('position %s out of the reference %s'%(position,total_length))
    return [seq, SNP_output]

def modelSNP(seq,Chr,num_mut_chr,num_indel_chr):
    total_length = len(seq)
    # indel modelling
    indel_output = []
    SNP_output = []
    seq = list(seq)
    if num_indel_chr > 0:
        candidate_position = [i for i in range(0, total_length) if seq[i] not in ['-','N']]
        indel_set = random.sample(candidate_position, k=num_indel_chr)
        seq, indel_output = modelindel(seq,Chr,indel_set)
    if num_mut_chr > 0:
        total_length = len(seq)
        # SNP modelling
        candidate_position = [i for i in range(0, total_length) if seq[i] not in ['-', 'N']]
        num_mut_chr = min(num_mut_chr,len(candidate_position))
        position_set = random.sample(candidate_position, k=num_mut_chr)
        SNP_output = []
        for position in position_set:
            gene_info = contig_to_gene(Chr, position)
            REF = seq[position]
            temp_ALT = ['A', 'T', 'G', 'C']
            try:
                temp_ALT.remove(REF)
            except ValueError:
                pass
            ALT = random.choices(temp_ALT, k=1)[0]
            seq[position] = ALT
            temp_line = [Chr,str(position+1),REF,ALT,'Other','None']
            if False and gene_info != []:
                # a gene
                Chr_gene, position_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
                codon_start = position_gene - 1 - int((position_gene - 1) % 3)
                if codon_start <= position_gene - 1:
                    Ref_seq_chr = Ref_seq[Chr_gene]
                    SNP_seq_chr = Ref_seq_chr
                    Ref_seq_chr = causeSNP(Ref_seq_chr, position_gene, REF,Reverse_chr)
                    Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                    if len(Ref_seq_codon) == 3:
                        Ref_seq_aa = translate(Ref_seq_codon)[0]
                        SNP_seq_chr = causeSNP(SNP_seq_chr, position_gene, ALT, Reverse_chr)
                        SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_aa = translate(SNP_seq_codon)[0]
                        temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                        temp_line[-1]=''.join([Ref_seq_aa,SNP_seq_aa])
                        temp_line[-2]=temp_NorS
            SNP_output.append('\t'.join(temp_line)+'\n')
    return [''.join(seq),SNP_output,indel_output]

def bamconverg(tempbamoutput):
    cmds = ''
    if not time_evaluation and run_bcf:
        try:
            f1 = open(tempbamoutput + '.sorted.sam','r')
        except IOError:
            cmds += '/usr/bin/time -v samtools sort -O SAM -@ %s %s.sam -o %s.sorted.sam\n' % (
                threadstouse, tempbamoutput,
                tempbamoutput)
    return cmds

def merge_sample(database,vcfoutput,tempbamoutput):
    cmds = ''
    if not time_evaluation and run_bcf:
        cmds += 'conda activate bt2\n'
        allsam = tempbamoutput + '.sorted.sam'
        if 'bcfdefault' in vcfoutput:
            # -Ov vcf format
            cmds += '/usr/bin/time -v bcftools mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -Ou -f %s %s | bcftools call -Ov -c --ploidy 1 --threads %s > %s.raw.vcf\n' % (
                threadstouse, database,
                allsam, threadstouse, vcfoutput)
        elif 'indel' in vcfoutput:
            cmds += '/usr/bin/time -v bcftools mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -Ou -A --indels-2.0 -f %s %s | bcftools call -Ov -A -M -c --ploidy 1 --threads %s > %s.raw.vcf\n' % (
                threadstouse, database,
                allsam, threadstouse, vcfoutput)
        else:
            cmds += '/usr/bin/time -v bcftools mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -Ou -A -f %s %s | bcftools call -Ov -A -M -c --ploidy 1 --threads %s > %s.raw.vcf\n' % (
                threadstouse, database,
                allsam, threadstouse, vcfoutput)
        if not time_evaluation:
            cmds += '/usr/bin/time -v bcftools view -H -v snps -i \'QUAL>=20 && MIN(DP)>=3\' %s.raw.vcf > %s.flt.snp.vcf \n' % (
                    vcfoutput, vcfoutput)
            cmds += '/usr/bin/time -v bcftools view -H -v indels -i \'QUAL>=20 && MIN(DP)>=3\' %s.raw.vcf > %s.flt.indel.vcf \n' % (
                    vcfoutput, vcfoutput)
            if 'ref.fasta' not in database:
                # keep raw.vcf for BaFr
                cmds += 'nohup rm %s.raw.vcf &\n' % (vcfoutput)
            #cmds += 'nohup rm %s %s.sam &\n' % (allsam, tempbamoutput)
        if 'MGSW' in tempbamoutput or 'MGBIG' in tempbamoutput: # MG SARS-COV-2
            cmds = cmds.replace('bcftools mpileup', 'bcftools mpileup -d300000')
    if time_evaluation:
        cmds += '/usr/bin/time -v samtools sort -O SAM -@ %s %s.sam -o %s.sorted.sam\n' % (
            threadstouse, tempbamoutput,
            tempbamoutput)
        cmds += 'conda activate bt2\n'
        cmds += '/usr/bin/time -v bcftools mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -Ou -A -f %s %s.sorted.sam | bcftools call -Ov -A -M -c --ploidy 1 --threads %s > %s.raw.vcf\n' % (
            threadstouse, database,
            tempbamoutput, threadstouse, vcfoutput)
        cmds += '#bcftools view -H -v snps -i \'QUAL>=20 && MIN(DP)>=3\' %s.raw.vcf > %s.flt.snp.vcf \n' % (
            vcfoutput, vcfoutput)
        cmds += '#bcftools view -H -v indels -i \'QUAL>=20 && MIN(DP)>=3\' %s.raw.vcf > %s.flt.indel.vcf \n' % (
            vcfoutput, vcfoutput)
        cmds += '#mv %s.flt* %s/merge/\n'%(vcfoutput,output_dir)
        cmds += 'rm %s.raw.vcf %s.flt*\n' % (vcfoutput,vcfoutput)
        cmds += 'rm %s.*\n' % (tempbamoutput)
    return cmds

def cleanup(tempbamoutput):
    cmds = ''
    if not time_evaluation:
        cmds += 'nohup rm %s* &\n' % (tempbamoutput)
    if 'covid' in tempbamoutput:
        cmds += 'nohup rm %s* &\n' % (tempbamoutput + '.sam')
    return cmds

def time_evaluation_vcfer(tempbamoutput,database):
    cmds = ''
    newoutput = '%s.mappersamtovcf.vcf' % (tempbamoutput.replace('/bwa/', '/merge/'))
    cmds += 'cp %s.sam %s.sam\n' % (tempbamoutput.replace(output_dir_temp, output_dir), tempbamoutput)
    cmds += '/usr/bin/time -v samtools sort -n -O SAM -@ %s %s.sam -o %s.sorted.sam\n' % (
        threadstouse, tempbamoutput,
        tempbamoutput)
    cmds += '/usr/bin/time -v java -Xms40g -Xmx40g -jar %s --num-threads %s --reference %s --in-sam %s.sorted.sam --out-vcf %s\n' % (
        latest_mappersamtovcf, threadstouse, database, tempbamoutput, newoutput)
    cmds += '#mv %s %s/merge/\n' % (newoutput, output_dir)
    cmds += 'rm -f %s.sorted.sam %s\n' % (tempbamoutput, newoutput)
    return cmds

def time_evaluation_parallel(tempbamoutput,database):
    cmds = ''
    newoutput = '%s.mappersamtovcf.vcf' % (tempbamoutput.replace('/bwa/', '/merge/'))
    cmds += 'cp %s.sam %s.sam\n' % (tempbamoutput.replace(output_dir_temp, output_dir), tempbamoutput)
    cmds += '/usr/bin/time -v samtools sort -n -O SAM -@ %s %s.sam -o %s.sorted.sam\n' % (
        threadstouse, tempbamoutput,
        tempbamoutput)
    for threadstousesub in list(range(1,11)) + [15,20,25,30]:
        cmds += '/usr/bin/time -v java -Xms40g -Xmx40g -jar %s --num-threads %s --reference %s --in-sam %s.sorted.sam --out-vcf %s\n' % (
            latest_mappersamtovcf, threadstousesub, database, tempbamoutput, newoutput)
        cmds += '/usr/bin/time -v java -Xms40g -Xmx40g -jar %s --num-threads %s --reference %s --in-sam %s.sorted.sam --out-vcf %s\n' % (
            latest_mappersamtovcfold, threadstousesub, database, tempbamoutput, newoutput)
    cmds += '#mv %s %s/merge/\n' % (newoutput, output_dir)
    cmds += 'rm -f %s.sorted.sam %s\n' % (tempbamoutput, newoutput)
    return cmds
def run_vcfer(tempbamoutput,files,files2,database):
    cmds = ''
    newoutput = '%s.mappersamtovcf.vcf' % (tempbamoutput.replace('/bwa/', '/merge/'))
    try:
        f1 = open(newoutput, 'r')
    except IOError:
        try:
            f1 = open('%s.sorted.sam' % (tempbamoutput), 'r')
        except IOError:
            cmds += '/usr/bin/time -v samtools sort -n -O SAM -@ %s %s.sam -o %s.sorted.sam\n' % (
                threadstouse, tempbamoutput,
                tempbamoutput)
        cmds += '/usr/bin/time -v java -Xms40g -Xmx40g -jar %s --num-threads %s --reference %s --in-sam %s.sorted.sam --out-vcf %s\n' % (
            latest_mappersamtovcf, threadstouse, database, tempbamoutput, newoutput)
        cmds += 'nohup rm -f %s.sorted.sam &\n' % (tempbamoutput)
    return cmds
def run_bowtie(files,files2,database,tempbamoutput):
    # generate code
    cmds = ''
    if penalty_test:
        # not using --ma for end to end model, --ignore-quals for mismatch quality, always use the highest penalty
        cmds += 'conda activate bt\nbowtie2 --threads %s --mp %s,%s --rfg %s,%s --rdg %s,%s --np %s --ignore-quals -x %s -1 %s -2 %s -S %s.sam\n' % (
            threadstouse, penalty[0],penalty[0],penalty[1],penalty[2],penalty[1],penalty[2],penalty[4], database, files, files2, tempbamoutput)
    elif parallel_test:
        cmds += time_evaluation_parallel(tempbamoutput,database)
    elif time_evaluation:
       cmds += time_evaluation_vcfer(tempbamoutput,database)
    elif mappertovcfonly:
        cmds += run_vcfer(tempbamoutput, files, files2, database)
    else:
        try:
            f1 = open('%s.sam' % (tempbamoutput), 'r')
        except IOError:
            cmds += 'conda activate bt\nbowtie2-build %s %s\n' % (database, database)
            cmds += 'conda activate bt\n/usr/bin/time -v bowtie2 --threads %s -x %s -1 %s -2 %s -S %s.sam\n' % (
                threadstouse, database, files, files2, tempbamoutput)
            cmds += run_vcfer(tempbamoutput, files, files2, database)
    return cmds

def run_minimap(files,files2,database,tempbamoutput):
    cmds = ''
    if penalty_test:
        # -A match score set as 2, -B -= 2 and --score-N += 2
        cmds += 'conda activate bt\nminimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s %s.mmi %s %s >%s.sam\n' % (
            threadstouse, penalty[0] - 2,penalty[1],penalty[2],penalty[3] + 2 ,penalty[3]-penalty[4] + 2, database, files, files2, tempbamoutput)
    elif parallel_test:
        cmds += time_evaluation_parallel(tempbamoutput,database)
    elif time_evaluation:
        cmds += time_evaluation_vcfer(tempbamoutput, database)
    elif mappertovcfonly:
        cmds += run_vcfer(tempbamoutput, files, files2, database)
    else:
        try:
            f1 = open('%s.sam' % (tempbamoutput), 'r')
        except IOError:
            cmds += 'conda activate bt\nminimap2 -d %s.mmi %s \n' % (database, database)
            cmds += 'conda activate bt\n/usr/bin/time -v minimap2 -ax sr -t %s %s.mmi %s %s >%s.sam\n' % (
                threadstouse,
                database, files, files2, tempbamoutput)
            cmds += run_vcfer(tempbamoutput, files, files2, database)
    return cmds

def run_bwa(files,files2,database,tempbamoutput):
    cmds = ''
    if penalty_test: 
        # -A match score set as 2, -B -= 2
        cmds += 'conda activate bt\nbwa mem -t %s -B %s -O %s -E %s -A %s %s %s %s > %s.sam\n' % (
            threadstouse, penalty[0] - 2,penalty[1],penalty[2],penalty[3]+ 2, database, files, files2, tempbamoutput)
    elif parallel_test:
        cmds += time_evaluation_parallel(tempbamoutput,database)
    elif time_evaluation:
        cmds += time_evaluation_vcfer(tempbamoutput, database)
    elif mappertovcfonly:
        cmds += run_vcfer(tempbamoutput, files, files2, database)
    else:
        try:
            f1 = open('%s.sam' % (tempbamoutput), 'r')
        except IOError:
            cmds += 'conda activate bt\nbwa index %s\n' % (database)
            cmds += 'conda activate bt\n/usr/bin/time -v bwa mem -t %s %s %s %s > %s.sam\n' % (
                threadstouse, database, files, files2, tempbamoutput)
            cmds += run_vcfer(tempbamoutput, files, files2, database)
    return cmds

def modelSNPall(Input_seq, Input_id, Length,num_mut,database_name,m):
    Output = []
    Output_SNP = []
    Output_indel = []
    chr_set_indel = random.choices(Input_id, k=indel_time, weights=Length)
    unique_chr_set_indel = list(set(chr_set_indel))
    total_seq = 0
    for chr in Input_id:
        # mutated proportional to length
        seq = Input_seq[chr]
        total_seq += len(seq)
    if num_mut < 1:
        total_num_mut = int(num_mut*total_seq) + 1
    else:
        total_num_mut = int(num_mut) + 1
    chr_set_snp = random.choices(Input_id, k=total_num_mut, weights=Length)
    unique_chr_set_snp = list(set(chr_set_snp))
    #print(total_seq, num_mut, total_num_mut,chr_set_snp)
    for chr in Input_id:
        seq = Input_seq[chr]
        num_indel_chr = 0
        if chr in unique_chr_set_indel:
            num_indel_chr = chr_set_indel.count(chr)
        num_snp_chr = 0
        if chr in unique_chr_set_snp:
            num_snp_chr = chr_set_snp.count(chr)
        newseq, newoutput, newoutputindel = modelSNP(seq, chr,num_snp_chr,num_indel_chr)
        Output_SNP += newoutput
        Output_indel += newoutputindel
        Output.append('>%s\n%s\n' % (chr,
                                     newseq))
    # output mutated genome
    output_fasta = os.path.join(output_dir, 'data/%s.%s.SNP%s.fasta' % (database_name, num_mut_name,m))
    f1 = open(output_fasta, 'w')
    f1.write(''.join(Output))
    f1.close()
    f1 = open(output_fasta + '.snp.txt', 'w')
    f1.write(''.join(Output_SNP))
    f1.close()
    f1 = open(output_fasta + '.indel.txt', 'w')
    f1.write(''.join(Output_indel))
    f1.close()
    print('done %s mutations in %s'%(num_mut,database_name))
    return output_fasta

# matching database to fastq
genome_fastq_mapping = dict()
if 'hybrid' in genome_root:
    allgenome = glob.glob('%s/*/*%s'%(genome_root,genome_name)) + glob.glob('%s/*/*/*%s'%(genome_root,genome_name))
    for genome in allgenome:
        genome_fastq_mapping.setdefault(genome,glob.glob('%s/*%s'%(os.path.split(genome)[0],fastq_name)))
# load database
if genome_fastq_mapping == {}:
    allfastq = glob.glob('%s/*%s'%(genome_root,fastq_name))
    print('%s/*%s' % (genome_root,genome_name))
    database = glob.glob('%s/*%s' % (genome_root,genome_name))[0]
    print(database)
    for fastq_file in allfastq:
        if 'ref.fasta' not in database:
            database = '%s/%s' % (genome_root, os.path.split(fastq_file)[-1].replace(fastq_name, genome_name))
            genome_fastq_mapping.setdefault(database, [fastq_file])
        else:
            genome_fastq_mapping.setdefault(database, [])
            genome_fastq_mapping[database].append(fastq_file)

print(genome_fastq_mapping)
for database in genome_fastq_mapping:
    database_name = os.path.split(database)[-1]
    # cause SNP
    if mutation_repeat:
        mutation_time = [str(x) for x in range(0,5)]
    else:
        mutation_time = ['']
    # Set up mutation rate
    mut_set = [0, 1e-6,2e-6,3e-6,4e-6,5e-6,1e-5,2e-5,3e-5,4e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,3e-2,4e-2,5e-2] #mutation rate [0, 5, 10, 50, 100, 1000, 5000, 10000, 50000, 100000, 200000]
    if parallel_test:
        mut_set = [1e-2]
    indel_time = 200
    indel_nonorf = list(range(2, 20))
    indel_nonorf.extend(list(range(-20, -2)))
    if 'SRR10971381' in database:
        # covid
        mut_set = [5, 10, 50, 100, 150, 200, 250, 300]
        indel_time = 20
        indel_nonorf = list(range(2, 20))
        indel_nonorf.extend(list(range(-20, -2)))
    if 'ref.fasta' in database:
        # simulation and real data BaFR
        mut_set = [0]
    indel_orf = indel_nonorf
    if correct_genome:
        mut_set = [0]
    print(database, mut_set, indel_time)
    for m in mutation_time:
        if len(mut_set) != 0:
            mut_time = len(mut_set)
        while mut_time > 0:
            num_mut = mut_set[mut_time - 1]
            num_mut_name = num_mut
            if 'covid' not in output_dir:
                num_mut_name = '%.0e'%(num_mut)
            # cause SNP
            if num_mut == 0:
                # corrected genome and control
                mutated_genome = os.path.join(output_dir, 'data/%s.%s.SNP%s.fasta' % (database_name, 0,m))
                try:
                    ftry = open(mutated_genome,'r')
                except IOError:
                    os.system('cp %s %s'%(database,mutated_genome))
                    f1 = open(mutated_genome + '.snp.txt', 'w')
                    f1.close()
            elif cause_SNP:
                # simulate mutated strain genomes
                database_file = (database.replace(genome_name,'.fasta.fna'))
                try:
                    open(database_file, 'r')
                except IOError:
                    os.system('prodigal -q -i %s -d %s' % (database, database_file))
                try:
                    mutated_genome = os.path.join(output_dir, 'data/%s.%s.SNP%s.fasta' % (database_name, num_mut_name, m))
                    ftry = open(mutated_genome,'r')
                except IOError:
                    Ref_seq, Length, Input_seq, Input_id = loaddatabase(database_file, database)
                    mutated_genome = modelSNPall(Input_seq, Input_id, Length,num_mut,database_name,m)
            else:
                mutated_genome = os.path.join(output_dir, 'data/%s.%s.SNP%s.fasta' % (database_name,num_mut_name, m))
            mutated_genome_filenamebase = os.path.split(mutated_genome)[-1]
            if mapping_file:
                for fastq_file in genome_fastq_mapping[database]:
                    # find fastq2
                    fastq_file2 = fastq_file.replace(fastq_name, fastq_name2)
                    mutated_genome_filename = mutated_genome_filenamebase
                    fastq_filename = os.path.basename(fastq_file).split(fastq_name)[0]
                    if 'ref.fasta' in database:
                        mutated_genome_filename = fastq_filename + '__' + mutated_genome_filename
                    # call SNPs by time bowtie2
                    cmds = ''
                    cmds += run_bowtie(fastq_file, fastq_file2,
                                          mutated_genome,
                                          os.path.join(output_dir_temp + '/bwa',
                                                       mutated_genome_filename + '.bowtie'))
                    f1 = open(os.path.join(input_script_sub, '%s.bowtiemapper.vcf.sh' % (mutated_genome_filename)), 'w')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                    f1.close()
                    cmds = ''
                    outputvcf = os.path.join(output_dir_temp + '/merge',
                                             mutated_genome_filename + '.bowtie')
                    samoutput = os.path.join(output_dir_temp + '/bwa',
                                             mutated_genome_filename + '.bowtie')
                    if not time_evaluation:
                        try:
                            f1 = open(outputvcf + '.indel.flt.snp.vcf', 'r')
                        except IOError:
                            cmds += merge_sample(mutated_genome, outputvcf + '.indel', samoutput)
                        try:
                            f1 = open(outputvcf + '.bcfdefault.flt.snp.vcf', 'r')
                        except IOError:
                            cmds += merge_sample(mutated_genome, outputvcf + '.bcfdefault',
                                                 samoutput)
                    try:
                        f1 = open(outputvcf + '.flt.snp.vcf','r')
                    except IOError:
                        cmds += merge_sample(mutated_genome, outputvcf, samoutput)
                    cmdssam = ''
                    if cmds != '':
                        cmdssam = bamconverg(samoutput)
                    cmds += cleanup(samoutput)
                    f1 = open(os.path.join(input_script_sub, '%s.bowtie.vcf.sh' % (mutated_genome_filename)), 'w')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmdssam + cmds))
                    f1.close()
                    # call SNPs by time bwa
                    cmds = ''
                    cmds += run_bwa(fastq_file, fastq_file2,
                                          mutated_genome,
                                          os.path.join(output_dir_temp + '/bwa',
                                                       mutated_genome_filename + '.bwa'))
                    f1 = open(os.path.join(input_script_sub, '%s.bwamapper.vcf.sh' % (mutated_genome_filename)), 'w')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                    f1.close()
                    cmds = ''
                    outputvcf = os.path.join(output_dir_temp + '/merge',
                                             mutated_genome_filename + '.bwa')
                    samoutput = os.path.join(output_dir_temp + '/bwa',
                                                              mutated_genome_filename + '.bwa')
                    if not time_evaluation:
                        try:
                            f1 = open(outputvcf + '.indel.flt.snp.vcf', 'r')
                        except IOError:
                            cmds += merge_sample(mutated_genome, outputvcf + '.indel',
                                                 samoutput)
                        try:
                            f1 = open(outputvcf + '.bcfdefault.flt.snp.vcf', 'r')
                        except IOError:
                            cmds += merge_sample(mutated_genome, outputvcf + '.bcfdefault',
                                                 samoutput)
                    try:
                        f1 = open(outputvcf + '.flt.snp.vcf', 'r')
                    except IOError:
                        cmds += merge_sample(mutated_genome, outputvcf, samoutput)
                    cmdssam = ''
                    if cmds != '':
                        cmdssam = bamconverg(samoutput)
                    cmds += cleanup( samoutput)
                    f1 = open(os.path.join(input_script_sub, '%s.bwa.vcf.sh' % (mutated_genome_filename)), 'w')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmdssam + cmds))
                    f1.close()
                    # call SNPs by time minimap
                    cmds = ''
                    cmds += run_minimap(fastq_file, fastq_file2,
                                          mutated_genome,
                                          os.path.join(output_dir_temp + '/bwa',
                                                       mutated_genome_filename + '.minimap'))
                    f1 = open(os.path.join(input_script_sub, '%s.minimapmapper.vcf.sh' % (mutated_genome_filename)), 'w')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                    f1.close()
                    cmds = ''
                    outputvcf = os.path.join(output_dir_temp + '/merge',
                                             mutated_genome_filename + '.minimap')
                    samoutput = os.path.join(output_dir_temp + '/bwa',
                                             mutated_genome_filename + '.minimap')
                    if not time_evaluation:
                        try:
                            f1 = open(outputvcf + '.indel.flt.snp.vcf', 'r')
                        except IOError:
                            cmds += merge_sample(mutated_genome, outputvcf + '.indel',
                                                 samoutput)
                        try:
                            f1 = open(outputvcf + '.bcfdefault.flt.snp.vcf', 'r')
                        except IOError:
                            cmds += merge_sample(mutated_genome, outputvcf + '.bcfdefault',
                                                 samoutput)
                    try:
                        f1 = open(outputvcf + '.flt.snp.vcf', 'r')
                    except IOError:
                        cmds += merge_sample(mutated_genome, outputvcf, samoutput)
                    cmdssam = ''
                    if cmds != '':
                        cmdssam = bamconverg(samoutput)
                    cmds += cleanup(samoutput)
                    f1 = open(os.path.join(input_script_sub, '%s.minimap.vcf.sh' % (mutated_genome_filename)), 'w')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmdssam + cmds))
                    f1.close()
                    # run bowtie and mapper on the same node
                    cmds = ''
                    cmds += 'bash %s\n' % (
                        os.path.join(input_script_sub, '%s.bowtiemapper.vcf.sh' % (mutated_genome_filename)))
                    cmds += 'bash %s\n' % (
                        os.path.join(input_script_sub, '%s.minimapmapper.vcf.sh' % (mutated_genome_filename)))
                    cmds += 'bash %s\n' % (
                        os.path.join(input_script_sub, '%s.bwamapper.vcf.sh' % (mutated_genome_filename)))
                    cmds += 'bash %s\n'%(os.path.join(input_script_sub, '%s.bowtie.vcf.sh' % (mutated_genome_filename)))
                    cmds += 'bash %s\n' % (os.path.join(input_script_sub, '%s.minimap.vcf.sh' % (mutated_genome_filename)))
                    cmds += 'bash %s\n' % (os.path.join(input_script_sub, '%s.bwa.vcf.sh' % (mutated_genome_filename)))
                    f1 = open(os.path.join(input_script_sub, '%s.bcftools.vcf.sh' % (mutated_genome_filename)), 'w')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                    f1.close()
            mut_time -= 1

if threadstouse < 10 and not time_evaluation:
    f1 = open(os.path.join(input_script, 'allsnpmodel1.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    f2 = open(os.path.join(input_script, 'allsnpmodel2.sh'), 'w')
    f2.write('#!/bin/bash\nsource ~/.bashrc\n')
    f3 = open(os.path.join(input_script, 'allsnpmodel3.sh'), 'w')
    f3.write('#!/bin/bash\nsource ~/.bashrc\n')
    f4 = open(os.path.join(input_script, 'allsnpmodel4.sh'), 'w')
    f4.write('#!/bin/bash\nsource ~/.bashrc\n')
    f5 = open(os.path.join(input_script, 'allsnpmodel5.sh'), 'w')
    f5.write('#!/bin/bash\nsource ~/.bashrc\n')
    total_test = 1
    if time_evaluation:
        total_test = 1 # run each pipeline 10 times
    i = 0
    for m in range(0,total_test):
        for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.bcftools.vcf.sh')):
            os.system('cp %s %s%s'%(sub_scripts,sub_scripts,m))
            i += 1
            if i % 5 == 1:
                f1.write('bash %s%s 2> %s%s.err 1> %s%s.out\n' % (sub_scripts, m, sub_scripts, m, sub_scripts, m))
            elif i % 5 == 2:
                f2.write('bash %s%s 2> %s%s.err 1> %s%s.out\n' % (sub_scripts, m, sub_scripts, m, sub_scripts, m))
            elif i % 5 == 3:
                f3.write('bash %s%s 2> %s%s.err 1> %s%s.out\n' % (sub_scripts, m, sub_scripts, m, sub_scripts, m))
            elif i % 5 == 4:
                f4.write('bash %s%s 2> %s%s.err 1> %s%s.out\n' % (sub_scripts, m, sub_scripts, m, sub_scripts, m))
            else:
                f5.write('bash %s%s 2> %s%s.err 1> %s%s.out\n' % (sub_scripts, m, sub_scripts, m, sub_scripts, m))

    os.system('rm %s/*.bcftools.vcf.sh'%(input_script_sub))
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
else:
    m = 0
    f1 = open(os.path.join(input_script, 'allsnpmodelall.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.bcftools.vcf.sh')):
        os.system('cp %s %s%s' % (sub_scripts, sub_scripts, m))
        f1.write('bash %s%s 2> %s%s.err 1> %s%s.out\n' % (sub_scripts, m, sub_scripts, m, sub_scripts, m))
    os.system('rm %s/*.bcftools.vcf.sh' % (input_script_sub))
    f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
