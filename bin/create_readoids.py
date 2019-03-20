#! /usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
import subprocess


parser = OptionParser(usage="usage: %prog [options] <fasta file>")

parser.add_option('--length', '-l', action="store", default="100", type='int',dest="rd_ln",
                  help='Readoid length (default=%default)')
parser.add_option('--distance', '-d', action="store", default="100000", type='int',dest="distance",
                  help='Distance between readoids - insert size (default=%default)')
parser.add_option('--output', '-o', action="store", default=None, type='string', dest="output",
                  help='Output fasta file (default=%default)')
parser.add_option('--window_size', '-w', action="store", default="10000", type='int', dest="window_size",
                  help='Window size to slide across fasta file (default=%default)')

parser.add_option('--coverage','-c',action="store",type='int',dest="coverage",
                  help="Target Coverage Level")
parser.add_option('--quality','-q',action="store", type='int', dest="quality_level",
                  help="Add in silico errors 1 per N bases")


(options, args) = parser.parse_args()

# check inputs or die
if len(args) !=1:
    parser.error("Must supply fasta file.")
    sys.exit(-1)

def makeReadoids(input,output):
    try: 				# open up output
        out_fh = open(output, 'w')
        in_fh  = open(input, 'r')
        pair_count = 0
        for record in SeqIO.parse(in_fh,'fasta'):				# we are telling module it is fasta
            seq = record.seq
            id  = record.id
            ln  = len(seq)
            last_position = ln - options.distance
            rd_ln = options.rd_ln
            distance = options.distance
            window_size = options.window_size

            for i in range(0,last_position,window_size):
                read1     = seq[i:i + rd_ln]
                rd2_start = i + (distance - rd_ln)
                rd2_end   = rd2_start + rd_ln
                read2     = seq[rd2_start:rd2_start + rd_ln]
                read2_rc  = read2.reverse_complement()
                header1   = 'p' + str(pair_count) + '_' + id + '_' + str(i)	 # cant combine strings and ints, hence use of str()
                header2   = 'p' + str(pair_count) + '_' + id + '_' + str(rd2_end) # cant combine strings and ints, hence use of str()
                out_fh.write(">%s\n%s\n" % (header1, read1))
                out_fh.write(">%s\n%s\n" % (header2, read2_rc))
                pair_count += 1

    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1

    # see ~/dev/ScaffoldAssessment.pl makeReadoids perl sub
def convertQVtoProbability(QV):
    #return 10**((-QV)/float(10))
    return 100

def covertErrors2QV(errors, bases):
    return -10 * Math.log10(errors/bases)


def addSequencingErrors(sequence, prob):
    mutation=0
    bases = ["A","C","G","T"]

    #check each base
    for i in range(0, len(sequence)):
        keep_base = random.randint(0, prob)
        #mutate positions
        if not keep_base:
            #print sequence
            new_base = bases[random.randint(0,3)]
            while new_base == sequence[i]:
                new_base = bases[random.randint(0,3)]
            sequence = sequence.tomutable()
            sequence[i] = new_base
            sequence = sequence.toseq()
            mutation += 1

    return sequence, mutation

def makeReadoidsByCoverage(input,output,coverage):
    out_1_fh = open(output + ".1.fasta", 'w')
    out_2_fh = open(output + ".2.fasta", 'w')
    in_fh  = open(input, 'r')
    map = open(output +".map", 'w')
    insert_size = options.distance
    read_length = options.rd_ln
    total_bases = 0
    total_mutation = 0
    for entry in SeqIO.parse(in_fh,'fasta'):
        entry_seq = entry.seq
        id = entry.id
        entry_length = len(entry_seq)
        reads_needed = (coverage * entry_length) / read_length
        for i in range(0, reads_needed/2):

            #pick pairs at random start sites
            start_pos = random.randint(0,entry_length-(insert_size+1))
            mate1 = entry_seq[start_pos:start_pos + read_length]
            mate_start = start_pos + (insert_size - read_length)
            mate2 = entry_seq[mate_start:mate_start + read_length]

            #Add some errors to the reads
            if options.quality_level:
                total_bases += read_length*2
                mate1,mutation = addSequencingErrors(mate1, options.quality_level)
                total_mutation += mutation
                mate2,mutation = addSequencingErrors(mate2, options.quality_level)
                total_mutation += mutation

            #Only use one header so that resulting fastas can be converted to bam
            readoid_header   = id + '_' + str(i)

            # create map of read name to starting point
            map.write(readoid_header + "\t" + str(start_pos)+"\n")
            map.write(readoid_header + "\t" + str(mate_start)+"\n")

            # write reads to fasta file
            out_1_fh.write(">%s\n%s\n" % (readoid_header, mate1))
            out_2_fh.write(">%s\n%s\n" % (readoid_header, mate2.reverse_complement()))

    if options.quality_level:
        print"Target Error Rate: 1/" + str(options.quality_level) + " (Q"+ str(options.quality_level) +")"
        print "Actual Error Rate: 1/" + str(total_bases/total_mutation)

def command(cmd):
    print cmd
    #subprocess.check_call(cmd, shell=True)
    output = subprocess.check_output(cmd, shell=True)
    return output

def main():

    output_fasta=options.output
    input_fasta = args[0]
    if not output_fasta:
        output_fasta=re.sub(".fasta",".readoid.fasta",input_fasta)

    if options.coverage:
        makeReadoidsByCoverage(input_fasta,output_fasta,options.coverage)
        # convert fastas into bam
        command("read_format_converter.py -d FR -o "+output_fasta+".bam "+output_fasta+".1.fasta "+output_fasta+".2.fasta")
        command("shortreadtool.pl -B "+output_fasta+".bam -c -p -Q localhost -r "+input_fasta+" -w 1 -a")
    else:
        makeReadoids(input_fasta,output_fasta)


if __name__ == "__main__":			# this will call main and start the script
    sys.exit(main())

