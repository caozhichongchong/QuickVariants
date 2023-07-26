import glob
import os
import statistics
from datetime import datetime
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all results (ref data and vcfs)",
                      type=str, default='.',
                      metavar='.')
# optional input genome
optional.add_argument("-cluster",
                      help="a reference genome to run, default is all genomes",
                      type=str, default='',
                      metavar='genome1')
################################################## Definition ########################################################
args = parser.parse_args()
# set up cutoff
min_maf_for_call = .8 #Remove sample*candidate
min_cov = 6 # at least 6 reads mapped to POS
# mapper using only middle depth
middle_depth_only = True
total_indel = 200
################################################### Function ########################################################
def load_indel(bowtievcf,allindel,blank_indel,set_blank_indel = False):
    ref_indel = [] #unique ref indels that were detected
    Qualified_indel = []
    indel_result = [0,0,0]# TP, same POS diff indel, diff POS
    indel_result_set = []
    try:
        for lines in open(bowtievcf, 'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, USELESS, REF, INDELset = lines_set[:5]
            INDELset = INDELset.split(',')
            Depthset = [float(i) for i in lines_set[9].split(':')[-1].split(',')]
            totaldepth = sum(Depthset)
            if totaldepth >= min_cov:
                for i in range(1, len(Depthset)):
                    subdepth = Depthset[i]
                    INDEL = INDELset[i - 1]
                    CHRPOS = '%s\t%s' % (CHR, POS)
                    if subdepth >= min_cov and subdepth / totaldepth >= min_maf_for_call:
                        if CHRPOS not in blank_indel:
                            # skip blank indel
                            indel_result_set.append(CHRPOS)
                            # qualified indel
                            Qualified_indel, indel_result, ref_indel = compareindelallsimple(bowtievcf, allindel,
                                                                                             [CHR, POS, REF, INDEL,
                                                                                              subdepth, totaldepth],
                                                                                             Qualified_indel,
                                                                                             indel_result, ref_indel)
    except UnicodeDecodeError:
        for lines in open(bowtievcf, 'r',encoding='latin'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, USELESS, REF, INDELset = lines_set[:5]
            INDELset = INDELset.split(',')
            Depthset = [float(i) for i in lines_set[9].split(':')[-1].split(',')]
            totaldepth = sum(Depthset)
            if totaldepth >= min_cov:
                for i in range(1, len(Depthset)):
                    subdepth = Depthset[i]
                    INDEL = INDELset[i - 1]
                    CHRPOS = '%s\t%s' % (CHR, POS)
                    if subdepth >= min_cov and subdepth / totaldepth >= min_maf_for_call:
                        if CHRPOS not in blank_indel:
                            # skip blank indel
                            indel_result_set.append(CHRPOS)
                            # qualified indel
                            Qualified_indel, indel_result, ref_indel = compareindelallsimple(bowtievcf, allindel,
                                                                                             [CHR, POS, REF, INDEL,
                                                                                              subdepth, totaldepth],
                                                                                             Qualified_indel,
                                                                                             indel_result, ref_indel)
    f1 = open(bowtievcf + '.filtered','w')
    f1.write('CHR\tPOS\tREF\tINDEL\tDepth\tTotal_depth\tCompare\n')
    f1.write(''.join(Qualified_indel))
    f1.close()
    if set_blank_indel:
        blank_indel = indel_result_set
    return [blank_indel,[total_indel - len(ref_indel)] + indel_result] # FN = [total_indel - len(ref_indel)]

def compareindelsub(refstring, indelstring):
    matchingpair = 0
    indellen = len(indelstring)
    reflen = len(refstring)
    for i in range(1,indellen):
        if indelstring[(i-1):(i+1)] in refstring:
            matchingpair += 1
    if matchingpair >= min(reflen,indellen) -1 -1:
        return True
    return False

def compareindel(indelset1,indelset2):
    REF1, INDEL1 = indelset1
    REF2, INDEL2 = indelset2
    if '-' not in REF1:
        # deletion
        if len(INDEL2) < len(REF2):
            #if REF2 in REF1 or REF1 in REF2 or (len(REF1)>5 and (REF1[3:] in REF2 or REF1[:-3] in REF2)):
            if REF2 in REF1 or REF1 in REF2 or compareindelsub(REF1, REF2):
                # same deletion or a subset
                return True
        return False
    else:
        # insertion
        if len(INDEL2) > len(REF2):
            if INDEL2 in INDEL1 or INDEL1 in INDEL2 or compareindelsub(INDEL1, INDEL2):
            #if INDEL2 in INDEL1 or INDEL1 in INDEL2 or (len(INDEL1)>5 and (INDEL1[3:] in INDEL2 or INDEL1[:-3] in INDEL2)):
                # same insertion or a subset
                return True
        return False

# print(compareindel(['AGCCAAATACA','A'],['AGCCAA','']))
# print(compareindel(['----------','ACTGAATAAC'],['AG','AGACTGAATAAC']))
# print(compareindel(['----','CCGG'],['','GCCG']))
# print(compareindel(['C','CCAAGGCTAAG'],['CAAGGCTC','']))
# print(compareindel(['-------------','GGTTCGGACGTGG'],['','GGTTCGGACGTGG']))
# print(compareindel(['GAAAA','GAAA'],['AA','A']))
def load_refindel(refindel):
    allindel = dict()
    allindel_len = dict()
    for lines in open(refindel, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR, POS, INDEL,REF = lines_set
        POS = int(POS)
        allindel.setdefault(CHR,{})
        if '-' in REF:
            # NEED CHANGE LATER
            #POS = POS - len(REF) - 1
            indel_len =len(INDEL)
        else:
            indel_len = len(INDEL) - len(REF)
        allindel[CHR].setdefault(POS,[REF, INDEL])
        allindel_len[indel_len] = allindel_len.get(indel_len,0) + 1
    return [allindel,allindel_len]

def compareindelallsimple(bowtievcf,allindel,indelset,Qualified_indel,indel_result,ref_indel):
    CHR, POS,REF, INDEL,subdepth,totaldepth = indelset
    POS = int(POS)
    Output = False
    if CHR in allindel:
        allindelCHR = allindel[CHR]
        for POSREF in allindelCHR:
            REF_REF, INDEL_REF = allindelCHR[POSREF]
            if abs(POS - POSREF) <= 10:
                # indel position within 10 bp
                Output = True
                if compareindel([REF_REF, INDEL_REF], [REF, INDEL]):
                    indel_result[0] += 1
                    Qualified_indel.append(
                        '%s\t%s\t%s\t%s\t%s\t%s\tTP\n' % (CHR, POS, REF, INDEL, subdepth, totaldepth))
                    if [CHR,POSREF] not in ref_indel:
                        ref_indel.append([CHR,POSREF])
                else:
                    indel_result[1] += 1
                    Qualified_indel.append(
                        '%s\t%s\t%s\t%s\t%s\t%s\tFP_samePOS\n' % (
                            CHR, POS, REF, INDEL, subdepth, totaldepth))
    if not Output:
        indel_result[2] += 1
        Qualified_indel.append(
                '%s\t%s\t%s\t%s\t%s\t%s\tFP_diffPOS\n' % (CHR, POS, REF, INDEL, subdepth, totaldepth))
    indel_FN = []
    indel_FN2 = []
    indel_TP = []
    for CHR in allindel:
        allindelCHR = allindel[CHR]
        for POSREF in allindelCHR:
            if [CHR,POSREF] not in ref_indel:
                indel_FN.append('%s\t%s\t%s\n'%(CHR,POSREF,'\t'.join(allindelCHR[POSREF])))
                indel_FN2.append('%s\t%s\t\n' % (CHR, POSREF))
            else:
                indel_TP.append('%s\t%s\t%s\n' % (CHR, POSREF, '\t'.join(allindelCHR[POSREF])))
    f1 = open(bowtievcf + '.TP', 'w')
    f1.write(''.join(indel_TP))
    f1.close()
    f1 = open(bowtievcf+'.FN','w')
    f1.write(''.join(indel_FN))
    f1.close()
    f1 = open(bowtievcf + '.FN.temp', 'w')
    f1.write(''.join(indel_FN2))
    f1.close()
    os.system('grep -T -f %s %s --no-group-separator > %s' % (
        bowtievcf+'.FN.temp',
        bowtievcf,
        bowtievcf + '.FN.vcf'))
    # os.system('sort -k3 -n %s | sort -k2 > %s' %
    #           (bowtievcf + '.temp', bowtievcf + '.FN.vcf' )
    #           )
    os.system('rm -rf %s %s' % (bowtievcf + '.temp',bowtievcf+'.FN.temp'))
    return [Qualified_indel,indel_result,ref_indel]

def compareindelall(allindel,indelset,Qualified_indel,indel_result, indel_result_set,ref_indel,blank_indel,indel_result_set2):
    CHR, POS,REF, INDEL,subdepth,totaldepth = indelset
    CHRPOS = '%s\t%s' % (CHR, POS)
    if CHRPOS not in blank_indel:
        # skip blank
        POS = int(POS)
        Output = False
        indel_len = len(INDEL) - len(REF)
        if CHR in allindel:
            allindelCHR = allindel[CHR]
            for POSREF in allindelCHR:
                if abs(POS - POSREF) <= 10:
                    # indel position within 10 bp
                    Output = True
                    REF_REF, INDEL_REF = allindelCHR[POSREF]
                    if compareindel([REF_REF, INDEL_REF], [REF, INDEL]):
                        if '-' in REF_REF:
                            indel_lenref = len(INDEL_REF)
                        else:
                            indel_lenref = len(INDEL_REF) - len(REF_REF)
                        indel_len = indel_lenref
                        indel_result[0] += 1
                        Qualified_indel.append(
                            '%s\t%s\t%s\t%s\t%s\t%s\tTP\n' % (CHR, POS, REF, INDEL, subdepth, totaldepth))
                        indel_result_set.setdefault(indel_len, [0, 0, 0])
                        indel_result_set[indel_len][0] += 1
                        if [CHR,POSREF] not in ref_indel:
                            ref_indel.append([CHR,POSREF])
                    else:
                        indel_result[1] += 1
                        Qualified_indel.append(
                            '%s\t%s\t%s\t%s\t%s\t%s\tFP_samePOS\n' % (
                                CHR, POS, REF, INDEL, subdepth, totaldepth))
                        indel_result_set.setdefault(indel_len, [0, 0, 0])
                        indel_result_set[indel_len][1] += 1
        if not Output:
            indel_result[2] += 1
            Qualified_indel.append(
                    '%s\t%s\t%s\t%s\t%s\t%s\tFP_diffPOS\n' % (CHR, POS, REF, INDEL, subdepth, totaldepth))
            indel_result_set.setdefault(indel_len, [0, 0, 0])
            indel_result_set[indel_len][2] += 1
        indel_result_set2.append(CHRPOS)
    else:
        print('indel %s in blank'%(CHRPOS))
    return [Qualified_indel,indel_result, indel_result_set, ref_indel,indel_result_set2]

def load_indelmapper(mappervcf, allindel,blank_indel,set_blank_indel = False):
    ref_indel = []
    Qualified_indel = []
    indel_result = [0, 0, 0]  # TP, same POS diff indel, diff POS
    indel_result_set = dict()
    insertion = ['', 0,'','',[],[]]#CHR POS REF ALT SUBdepth Totaldepth
    deletion = ['', [0],'','',[],[]]#CHR POS REF ALT SUBdepth Totaldepth
    indel_result_set2 = []
    for lines in open(mappervcf, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR, POS = lines_set[:2]
        if '-' in lines:
            POS = int(POS)
            INDELset = lines_set[3].split(',')
            # using only middle depth
            Depthsetend = []
            Depthset = [i for i in lines_set[5].split(';')]
            if not middle_depth_only:
                Depthsetend = [i for i in lines_set[6].split(';')]
                totaldepth = float(lines_set[4])
            else:
                totaldepth = sum([float(j.split(',')[0]) + float(j.split(',')[1]) for j in Depthset])
            REF = lines_set[2]
            if totaldepth >= min_cov:
                for i in range(1, len(Depthset)):
                    subdepth = sum([float(j) for j in Depthset[i].split(',')])
                    if not middle_depth_only:
                        subdepth += sum([float(j) for j in Depthsetend[i].split(',')])
                    INDEL = INDELset[i-1]
                    # extending indels, loose cutoff #need mentioning it in the paper
                    if subdepth >= min_cov and subdepth / totaldepth >= min_maf_for_call - 0.1:
                        # qualified indel
                        if POS < 0:
                            #print(insertion, REF, INDEL,subdepth, totaldepth)
                            POS = -POS
                            # insertion
                            if POS != insertion[1]:
                                # diff insertion
                                # sum the previous insertion
                                if len(insertion[3]) > 1:  # insertion >= 2bp
                                    insertion[4] = statistics.mean(insertion[4])  # mean of subdepth
                                    insertion[5] = statistics.mean(insertion[5])  # mean of depth
                                    Qualified_indel, indel_result, indel_result_set,ref_indel,indel_result_set2 = compareindelall(allindel,
                                                                                                      insertion,
                                                                                                      Qualified_indel,
                                                                                                      indel_result,
                                                                                                      indel_result_set,ref_indel,blank_indel,indel_result_set2)
                                # start of a new insertion, stringent cutoff
                                if subdepth / totaldepth >= min_maf_for_call:
                                    insertion = [CHR, POS, '', INDEL, [subdepth], [totaldepth]]
                                else:
                                    insertion = ['', 0, '', '', [], []]
                            else:
                                insertion[3] += INDEL
                                insertion[4].append(subdepth)
                                insertion[5].append(totaldepth)
                            #print(insertion)
                        elif INDEL == '-':
                            #print(deletion, CHR, POS, subdepth, totaldepth)
                            # deletion
                            if CHR != deletion[0] or abs(int(POS) - int(deletion[1][-1])) > 10:
                                # diff deletion
                                # sum the previous deletion
                                if len(deletion[2]) > 1: # deletion >= 2bp
                                    deletion[4] = statistics.mean(deletion[4]) # mean of subdepth
                                    deletion[5] = statistics.mean(deletion[5]) # mean of depth
                                    deletion[1] = deletion[1][0] # use the first POS
                                    Qualified_indel, indel_result,indel_result_set,ref_indel,indel_result_set2 = compareindelall(allindel, deletion,
                                                                                    Qualified_indel, indel_result, indel_result_set,ref_indel,blank_indel,indel_result_set2)
                                # start of a new deletion, stringent cutoff
                                if subdepth / totaldepth >= min_maf_for_call:
                                    deletion = [CHR,[POS], REF, '', [subdepth], [totaldepth]]
                                else:
                                    deletion = ['', [0], '', '', [], []]
                            else:
                                # same deletion
                                deletion[1].append(POS)
                                deletion[2] += REF
                                deletion[4].append(subdepth)
                                deletion[5].append(totaldepth)
                            #print(deletion)
    # sum the last deletion
    if len(deletion[2]) > 1:  # deletion >= 2bp
        deletion[4] = statistics.mean(deletion[4]) # mean of subdepth
        deletion[5] = statistics.mean(deletion[5]) # mean of depth
        deletion[1] = deletion[1][0]  # use the first POS
        Qualified_indel, indel_result, indel_result_set,ref_indel,indel_result_set2 = compareindelall(allindel, deletion,
                                                                          Qualified_indel, indel_result,
                                                                          indel_result_set,ref_indel,blank_indel,indel_result_set2)
    # sum the last insertion
    if len(insertion[3]) > 1:  # insertion >= 2bp
        insertion[4] = statistics.mean(insertion[4])  # mean of subdepth
        insertion[5] = statistics.mean(insertion[5])  # mean of depth
        Qualified_indel, indel_result, indel_result_set,ref_indel,indel_result_set2 = compareindelall(allindel, insertion,
                                                                          Qualified_indel, indel_result,
                                                                          indel_result_set,ref_indel,blank_indel,indel_result_set2)
    f1 = open(mappervcf + '.indel.vcf.filtered', 'w')
    f1.write('CHR\tPOS\tREF\tINDEL\tDepth\tTotal_depth\tCompare\n')
    f1.write(''.join(Qualified_indel))
    f1.close()
    for indel_len in indel_result_set:
        indel_len_result = indel_result_set[indel_len]
        f2.write('%s\t%s\tmapper\t%s\t%s\t%s\t%s\n' % (indel_len,samplename, indel_len_result[0], max(allindel_len.get(indel_len,0) - indel_len_result[0],0), indel_len_result[1], indel_len_result[2]))
    if set_blank_indel:
        blank_indel = indel_result_set2
    return [blank_indel,[total_indel - len(ref_indel)] + indel_result]

def indel_compare(allsum,samplename,blank_indel_set,set_blank_indel):
    blank_indelbowtie, blank_indelbwa, blank_indelminimap, \
        blank_indelbowtiebcfdefault, blank_indelbwabcfdefault, blank_indelminimapbcfdefault,blank_indelmapper,\
        blank_indelmappernoancestor, blank_indelmappernoxmer19, blank_indelmappernoxmer20, blank_indelmappernoxmer21,blank_indelmapperbcf,\
        blank_indelbowtiemappersamtovcf,blank_indelmappermappersamtovcf,blank_indelminimapmappersamtovcf,blank_indelbwamappersamtovcf = blank_indel_set
    # bowtie
    # bowtievcf = '%s/merge/%s.bowtie.flt.indel.vcf' % (resultdir, samplename)
    # print('%s process bowtie indel %s' % (datetime.now(), bowtievcf))
    # blank_indelbowtie, indel_result = load_indel(bowtievcf, allindel, blank_indelbowtie,set_blank_indel)
    # allsum.append('%s\tbowtie2\t%s\t%s\t%s\t%s\n' % (
    # samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # # bwa
    # bwavcf = '%s/merge/%s.bwa.flt.indel.vcf' % (resultdir, samplename)
    # print('%s process bwa indel %s' % (datetime.now(), bwavcf))
    # blank_indelbwa, indel_result = load_indel(bwavcf, allindel, blank_indelbwa,set_blank_indel)
    # allsum.append(
    #     '%s\tbwa\t%s\t%s\t%s\t%s\n' % (samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # # minimap
    # minimapvcf = '%s/merge/%s.minimap.flt.indel.vcf' % (resultdir, samplename)
    # print('%s process minimap indel %s' % (datetime.now(), minimapvcf))
    # blank_indelminimap, indel_result = load_indel(minimapvcf, allindel, blank_indelminimap,set_blank_indel)
    # allsum.append('%s\tminimap2\t%s\t%s\t%s\t%s\n' % (
    # samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # if 'covid' not in args.i:
    #     # bowtie
    #     bowtievcf = '%s/merge/%s.bowtie.bcfdefault.flt.indel.vcf' % (resultdir, samplename)
    #     print('%s process bowtie indel %s' % (datetime.now(), bowtievcf))
    #     blank_indelbowtiebcfdefault, indel_result = load_indel(bowtievcf, allindel, blank_indelbowtiebcfdefault, set_blank_indel)
    #     allsum.append('%s\tbowtie2 (bcfdefault)\t%s\t%s\t%s\t%s\n' % (
    #         samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    #     # bwa
    #     bwavcf = '%s/merge/%s.bwa.bcfdefault.flt.indel.vcf' % (resultdir, samplename)
    #     print('%s process bwa indel %s' % (datetime.now(), bwavcf))
    #     blank_indelbwabcfdefault, indel_result = load_indel(bwavcf, allindel, blank_indelbwabcfdefault, set_blank_indel)
    #     allsum.append(
    #         '%s\tbwa (bcfdefault)\t%s\t%s\t%s\t%s\n' % (samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    #     # minimap
    #     minimapvcf = '%s/merge/%s.minimap.bcfdefault.flt.indel.vcf' % (resultdir, samplename)
    #     print('%s process minimap indel %s' % (datetime.now(), minimapvcf))
    #     blank_indelminimapbcfdefault, indel_result = load_indel(minimapvcf, allindel, blank_indelminimapbcfdefault, set_blank_indel)
    #     allsum.append('%s\tminimap2 (bcfdefault)\t%s\t%s\t%s\t%s\n' % (
    #         samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # # mapper
    # mappervcf = '%s/merge/%s.mapper1.vcf.snp' % (resultdir, samplename)
    # print('%s process mapper indel %s' % (datetime.now(), mappervcf))
    # blank_indelmapper, indel_result = load_indelmapper(mappervcf, allindel, blank_indelmapper,set_blank_indel)
    # allsum.append('%s\tmapper\t%s\t%s\t%s\t%s\n' % (
    # samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # mapper no ancestor
    mappervcf = '%s/merge/%s.mapper1.noancestor.vcf.snp' % (resultdir, samplename)
    print('%s process mapper indel %s' % (datetime.now(), mappervcf))
    blank_indelmappernoancestor, indel_result = load_indelmapper(mappervcf, allindel, blank_indelmapper, set_blank_indel)
    allsum.append('%s\tmapper (noancestor)\t%s\t%s\t%s\t%s\n' % (
        samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # mapper no ancestor no x-mer
    mappervcf = '%s/merge/%s.mapper1.kmer15.vcf.snp' % (resultdir, samplename)
    print('%s process mapper indel %s' % (datetime.now(), mappervcf))
    blank_indelmappernoxmer19, indel_result = load_indelmapper(mappervcf, allindel, blank_indelmapper,
                                                                 set_blank_indel)
    allsum.append('%s\tmapper (15-mer)\t%s\t%s\t%s\t%s\n' % (
        samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # mapper no ancestor no x-mer
    mappervcf = '%s/merge/%s.mapper1.kmer20.vcf.snp' % (resultdir, samplename)
    print('%s process mapper indel %s' % (datetime.now(), mappervcf))
    blank_indelmappernoxmer20, indel_result = load_indelmapper(mappervcf, allindel, blank_indelmapper,
                                                             set_blank_indel)
    allsum.append('%s\tmapper (20-mer)\t%s\t%s\t%s\t%s\n' % (
        samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # mapper no ancestor no x-mer
    mappervcf = '%s/merge/%s.mapper1.kmer25.vcf.snp' % (resultdir, samplename)
    print('%s process mapper indel %s' % (datetime.now(), mappervcf))
    blank_indelmappernoxmer21, indel_result = load_indelmapper(mappervcf, allindel, blank_indelmapper,
                                                             set_blank_indel)
    allsum.append('%s\tmapper (25-mer)\t%s\t%s\t%s\t%s\n' % (
        samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # # mapper bcftools
    # mappervcf = '%s/merge/%s.mapper1.bcf.vcf.flt.indel.vcf' % (resultdir, samplename)
    # print('%s process mapper indel %s' % (datetime.now(), mappervcf))
    # blank_indelmapperbcf, indel_result = load_indel(mappervcf, allindel, blank_indelmapper,
    #                                                            set_blank_indel)
    # allsum.append('%s\tmapper (bcf)\t%s\t%s\t%s\t%s\n' % (
    #     samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # bowtie mapper sam to vcf
    bowtievcf = '%s/merge/%s.bowtie.mappersamtovcf.vcf.snp' % (resultdir, samplename)
    print('%s process bowtievcf indel %s' % (datetime.now(), bowtievcf))
    blank_indelbowtiemappersamtovcf, indel_result = load_indelmapper(bowtievcf, allindel, blank_indelmapper,
                                                             set_blank_indel)
    allsum.append('%s\tbowtie (mapper samtovcf)\t%s\t%s\t%s\t%s\n' % (
        samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # minimap mapper sam to vcf
    bowtievcf = '%s/merge/%s.minimap.mappersamtovcf.vcf.snp' % (resultdir, samplename)
    print('%s process minimap indel %s' % (datetime.now(), bowtievcf))
    blank_indelminimapmappersamtovcf, indel_result = load_indelmapper(bowtievcf, allindel, blank_indelmapper,
                                                                     set_blank_indel)
    allsum.append('%s\tminimap (mapper samtovcf)\t%s\t%s\t%s\t%s\n' % (
        samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    # bwa mapper sam to vcf
    bowtievcf = '%s/merge/%s.bwa.mappersamtovcf.vcf.snp' % (resultdir, samplename)
    print('%s process bwa indel %s' % (datetime.now(), bowtievcf))
    blank_indelbwamappersamtovcf, indel_result = load_indelmapper(bowtievcf, allindel, blank_indelmapper,
                                                                     set_blank_indel)
    allsum.append('%s\tbwa (mapper samtovcf)\t%s\t%s\t%s\t%s\n' % (
        samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    #mapper mapper sam to vcf
    bowtievcf = '%s/merge/%s.mapper1.mappersamtovcf.vcf.snp' % (resultdir, samplename)
    print('%s process mapper indel %s' % (datetime.now(), bowtievcf))
    blank_indelmappermappersamtovcf, indel_result = load_indelmapper(bowtievcf, allindel, blank_indelmapper,
                                                                     set_blank_indel)
    allsum.append('%s\tmapper (mapper samtovcf)\t%s\t%s\t%s\t%s\n' % (
        samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
    blank_indel_set = [blank_indelbowtie, blank_indelbwa, blank_indelminimap,
                       blank_indelbowtiebcfdefault, blank_indelbwabcfdefault, blank_indelminimapbcfdefault,blank_indelmapper,
                       blank_indelmappernoancestor, blank_indelmappernoxmer19,blank_indelmappernoxmer20,
                       blank_indelmappernoxmer21,blank_indelmapperbcf,
                       blank_indelbowtiemappersamtovcf,blank_indelmappermappersamtovcf,blank_indelminimapmappersamtovcf,blank_indelbwamappersamtovcf]
    return [allsum,blank_indel_set]

##################################################### Main ##########################################################
allresultdir = glob.glob('%s/SNP_model*/'%(args.i))
for resultdir in allresultdir:
    resultdirfilename = os.path.split(resultdir)[-1]
    allref = glob.glob('%s/data/%s*.5.*.indel.txt' % (resultdir, args.cluster)) + \
             glob.glob('%s/data/%s*.1e-01.*.indel.txt' % (resultdir, args.cluster))
    for ref in allref:
        try:
            f1 = open(ref.replace('.5.', '.0.').replace('.1e-01.', '.0.'), 'r')
        except IOError:
            f1 = open(ref.replace('.5.','.0.').replace('.1e-01.', '.0.'),'w')
            print('creating 0 SNP background genome',ref.replace('.5.','.0.').replace('.1e-01.', '.0.'))
            f1.close()
        try:
            f1 = open(ref.replace('.0.SNP.fasta.indel.txt', '.0.SNP.fasta.snp.txt'), 'r')
        except IOError:
            f1 = open(ref.replace('.0.SNP.fasta.indel.txt', '.0.SNP.fasta.snp.txt'), 'w')
            print('creating 0 SNP background genome',ref.replace('.5.','.0.'))
            f1.close()
    allref = glob.glob('%s/data/%s*corrected.0.*.indel.txt'%(resultdir,args.cluster))
    allsum = ['Sample\tTool\tFN\tTP\tFP_samePOS\tFP_diffPOS\n']
    f2 = open('%s/modelindelsumlen_%s.txt' % (resultdir,resultdirfilename), 'w')
    f2.write('Indellen\tSample\tTool\tTP\tFN\tFP_samePOS\tFP_diffPOS\n')
    # load blank indel at 0 SNPS -> not used because indel shifts
    Ref_set = dict()
    for ref in allref:
        print('%s load blank ref indel %s' % (datetime.now(), ref))
        allindel, allindel_len = load_refindel(ref)
        samplename = os.path.split(ref)[-1].split('.indel.txt')[0]
        genomename = samplename.split('.')[0]
        allsum, blank_indel_set = indel_compare(allsum,samplename,[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]],True)
        Ref_set.setdefault(genomename,blank_indel_set)
        print(genomename,[len(x) for x in blank_indel_set])
    # load all other indels
    allref = glob.glob('%s/data/%s*corrected*.indel.txt' % (resultdir, args.cluster))
    for ref in allref:
        if '.500000.' not in ref and '.0.' not in ref and '3e-05' not in ref and '3e-06' not in ref: # need change
            print('%s load ref indel %s' % (datetime.now(), ref))
            allindel, allindel_len = load_refindel(ref)
            samplename = os.path.split(ref)[-1].split('.indel.txt')[0]
            genomename = samplename.split('.')[0]
            blank_indel_set = Ref_set[genomename]
            print(genomename, [len(x) for x in blank_indel_set])
            allsum, blank_indel_set = indel_compare(allsum, samplename, blank_indel_set, False)
    f1 = open('%s/modelindelsum_%s.txt' % (resultdir,resultdirfilename), 'w')
    f1.write(''.join(allsum))
    f1.close()
    f2.close()
