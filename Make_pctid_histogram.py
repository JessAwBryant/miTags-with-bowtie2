__author__ = 'Jessica_Bryant'
"""
Input: a directory path containing already parsed bowtie2 output by John E's script filter_blast_m8.py
       all files need the extension 'miTags.db.L70.I.97.P'

Output: a histogram plotting the distribution of percent ids of recruited reads to database sequences

python /mnt/lysine/jbryant/DEEP_SEQUENCING/riboframe_bowtie/Make_pctid_histogram.py
"""

import glob as glob
from collections import defaultdict
import string
import matplotlib.pyplot as plt
import sys as sys
import re as re

def parseCigar(cigarString):
    #Daniels Function

    dictCigarCounts = defaultdict(int)
    listTypes = [i for i in cigarString if not i.isdigit()]
    inTable = "MIDSHPNX="
    outTable = "_"*len(inTable)
    translationTable = string.maketrans(inTable, outTable)
    listNumbers = cigarString.translate(translationTable).split("_")

    for position, type in enumerate(listTypes):
        count = listNumbers[position]
        dictCigarCounts[type] += int(count)

    seqLength = 0.0
    alignLength = 0.0
    alignCigarMismatches = 0.0
    alignLengthOnRef = 0.0

    for type in dictCigarCounts:
        if type in "MISH=X":
            seqLength += dictCigarCounts[type]
        if type in "MI=X":
            alignLength += dictCigarCounts[type]
        if type in "IDX":
            alignCigarMismatches += dictCigarCounts[type]
        if type in "MD=X":
            alignLengthOnRef += dictCigarCounts[type]
    a = alignLength
    pctd_id = 1 - (alignCigarMismatches/a)
    #print listTypes, listNumbers, seqLength, alignLength, alignCigarMismatches, alignLengthOnRef
    return(pctd_id)


if __name__ == '__main__':
    file_paths = sys.argv[1]

    excluded_runs = ['HOT233_1_0125m', 'HOT227_1_1000m', 'HOT233_1_0025m', 'HOT231_1_0125m', 'HOT227_1_0500m',
                     'HOT237_2_1000m', 'HOT231_1_0025m', 'HOT238_2_1000m', 'HOT234_2_0025m', 'HOT234_2_0075m',
                     'HOT237_1_0075m']

    l = []

    for file in glob.glob(file_paths+"/"+"*m.16S.fasta.bowtie.SILVA_123_miTags.db.L70.I.97.P"):
    #for file in glob.glob("HOT238_2_0500m.16S.fasta.bowtie.SILVA_123_miTags.db.L70.I.97.P"):
        if re.search('(HOT.*?m)', file).group(1) in excluded_runs:
            print 'skipped: ', re.search('(HOT.*?m)', file).group(1)
            continue
        print file
        file_handle=open(file)
        line = file_handle.readline()
        while line:
            cigar=line.split('\t')[5]
            l.append(parseCigar(cigar))
            #if float(parseCigar(cigar)) < 0.94: print cigar, parseCigar(cigar)
            if float(parseCigar(cigar)) < 0.94: print cigar, parseCigar(cigar)
            line = file_handle.readline()

    print 'total fragments: ', len(l)
    print 'pctds greater than or equal to 0.99: ', sum([1.0 for x in l if x >= 0.99])/len(l)

    # the histogram of the data
    n, bins, patches = plt.hist(l, 30, facecolor='green', alpha=0.75)

    plt.xlabel('identity to database sequence (%)')
    plt.ylabel('Counts')
    #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
    plt.axis([0.95, 1.0, 0, 3500000])
    plt.grid(True)

    plt.savefig("foo.pdf", bbox_inches='tight')


