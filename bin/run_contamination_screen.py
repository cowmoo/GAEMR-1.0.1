#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import subprocess
import sys
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()
from optparse import OptionParser
from gaemr.BlastFilter import VecScreenFilter, gContamFilter, MitoFilter, rRNAFilter, CombinedContaminationObject, BlastFilter
from gaemr.Md5 import Md5
from gaemr.Directive import Directive
from gaemr.DirectiveDetails import DirectiveDetails
from gaemr.PostProcessingFile import PostProcessingFile

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option('--input', '-i', action="store", default=None, type='string', dest="input",
    help='Input file fasta')

parser.add_option('--output', '-o', action="store", default="contam_screen", type='string', dest="output",
    help='Output file header: ')

parser.add_option('--output_format', '-f', action="store", default="m8", type='string', dest="output_format",
    help='Output file format: (default=%default)\n\tm8 - BLAST m8 format\n\tSubPrep - SubmissionPrep format contig_id, contig length, contam start, contam end, hit description\n\tgeneral - General removal format contig_name, contig length, contam start, contam end\n\tDirective - action, contig_id, contam start, contam_end, hit description')

parser.add_option('--stdout', '', action="store_true", default=False, dest="stdout",
    help='Output to stdout (default= False)')

parser.add_option('--redundant', '-r', action="store_true", default=False, dest="redundant",
    help='Remove redundant hits (default=%default)')

parser.add_option('--extend', '-e', action="store_true", default=False, dest="extend",
    help='Merge and extend overlapping hits (default=%default)')

parser.add_option('--reach', '-R', action="store", type='int', default=200, dest="reach",
    help='Size of reach (distance between hits) for extension and merging (default=%default)')

parser.add_option('--threads', '-t', action="store", default="4", type='int', dest="threads",
    help='Number of threads to use for blastn (default=%default)')

### All
parser.add_option('--all', '', action="store_true", default=False, dest="all",
    help='Run all screens: gcontam, mitochondria, common and Vecscreen (default=%default)')

### gContam
parser.add_option('--gcontam', '-g', action="store_true", default=False, dest="gcontam",
    help='Run Screen against NCBI\'s gcontam db (default=%default)')

### Mito
parser.add_option('--mito', '-m', action="store_true", default=False, dest="mito",
    help='Run mitochondria screen (default=%default)')

parser.add_option('--pct_mito', '-p', action="store", default=".30", type='float', dest="pct_mito",
    help='''Fraction of contig required to hit mitochondria to be filtered (default=0.30).
Does not apply if the --use_mito_len option is specified.''')

parser.add_option('--len_mito', '', action="store", default=120, type='int', dest="len_mito",
    help='''The lengh cuttoff to use for mitochondria screen.
Only applies with the --use_mito_len option is specified.''')

parser.add_option('--use_mito_len', '', action="store_true", default=False, dest="use_mito_len",
    help='Use a length cutoff for mitochondria screen.  Otherwise a percent coverage cutoff is used.')

### special processing if use a combined mito and gcontam db
parser.add_option('--no_gcontammitodb', '', action="store_false", default=True, dest="gcontammitodb",
    help='''Do not use the combined gcontam1 and mito DB for either the gcontam
or the mito screens.  The combined mitogcontam DB is used otherwise.''')

### VecScreen
parser.add_option('--vecscreen', '-v', action="store_true", default=False, dest="vecscreen",
    help='Run NCBI VecScreen (default=%default)')

### Common Contaminants (prok and euk)

parser.add_option('--common_prok', '-P', action="store_true", dest="common_prok", default=False,
                  help="Run common prokaryote contaminants screen.")

parser.add_option('--common_euk', '-E', action="store_true", dest="common_euk", default=False,
                  help="Run common eukaryote contaminants screen.")

### use NCBI special adaptors_for_screening db
parser.add_option('--adaptorsdb', '', action="store_true", default=False, dest="adaptorsdb",
    help='For VecScreen use apaptors_for_screening DB instead of UniVec DB which is the default')

parser.add_option('--no_protection', action="store_false", default=True, dest="use_protection",
    help='''For VecScreen do not protect weak hits using a default list of protected terms.
Protection is enabled by default''')

parser.add_option('--illumina_primer', '-I', action="store_true", default=False, dest="primer_only",
    help='Only report hits to Illumina adapter/primer in VecScreen (default=%default)')

### rRNA
parser.add_option('--rRNA', '', action="store_true", default=False, dest="rRNA",
    help='Run NCBI rRNA Screen (default=%default)')

###Custom
parser.add_option('--custom', '-c', action="store", default=None, type="string", dest="custom_db",
    help='Custom db for contamination screen')

parser.add_option('--filter_as', '-F', action="store", default=None, type="string", dest="filter_as",
    help='Filter hits from custom database like:\tgcontam\tmitochondria\tvecscreen\trRNA\tNone (default=%default)')

(options, args) = parser.parse_args()

# check inputs or die
if not options.input:
    parser.error("Must supply input file with --input.")
    sys.exit(-1)

if not options.output:
    parser.error("Must supply output file with --output.")
    sys.exit(-1)

if not options.mito and \
        not options.gcontam and \
        not options.vecscreen and \
        not options.all and \
        not options.custom_db and \
        not options.rRNA and \
        not options.common_prok and \
        not options.common_euk:
    parser.error("Must select a screen: --all, -g, -m, -v, -P, -E, --rRNA or -c")
    sys.exit(-1)

if options.all:
    options.mito = True
    options.gcontam = True
    options.vecscreen = True
    options.common_prok = True
    options.common_euk = True
    #options.rRNA = True

### Change pct_mito from percent to fraction
if options.pct_mito > 1:
    pct_mito = options.pct_mito / 100
else:
    pct_mito = options.pct_mito

def log_filter(filter, Name):
    hit_list = filter.get_removal_coords()
    print "Filter: " + Name + "\t" + str(len(hit_list))


def command(cmd):
    print cmd
    subprocess.check_call(cmd, shell=True)


def vecscreen_filters(xml):
    vec = VecScreenFilter(xml, use_protection=options.use_protection)
    vec.filter_weak_hits()
    if options.primer_only:
        vec.select_Illumina_hits()
    log_filter(vec, "Vector Screen")
    return vec


def rRNA_filters(xml):
    rRNA = rRNAFilter(xml)
    rRNA.filter_by_length(100)
    rRNA.filter_by_pct_covered(0.8)
    log_filter(rRNA, "rRNA Screen")
    return rRNA


def gcontam_filters(xml):
    gcontam = gContamFilter(xml)
    if options.gcontammitodb:
        gcontam.ignore_mito_hits()
    gcontam.filter_hits()
    log_filter(gcontam, "gContam Screen")
    return gcontam


def mito_filters(xml):
    mito = MitoFilter(xml, options.len_mito if options.use_mito_len else pct_mito)

    if options.use_mito_len:
        mito.select_hits_by_length()
    else:
        mito.select_valid_mito_hits()
    log_filter(mito, "Mitochondria Screen")
    return mito


def no_filter(xml):
    raw = BlastFilter(xml)
    log_filter(raw, "Raw Hits - No Filter")
    return raw


def list_to_string(list):
    out_string = ""
    for hit in list:
        out_string += (''.join([str(i) + "\t" for i in hit]))
        out_string = out_string[:-1] + "\n"

    return out_string


def __write_directives_file(combined):
        details_file = PostProcessingFile(options.output + ".details.txt")
        directives_file = PostProcessingFile(options.output + ".directives.txt")
        md5 = Md5()
        directives_file.add_output_line(Directive(action = "checksum", id = md5.get_md5_of_file_contents(options.input),
                                                  comment = "GAEMR Md5 checksum of " + options.input))
        
        for line_list in combined.get_removal_coords(False):
            directives_file.add_output_line(Directive(action = "remove", id = line_list[0], start = line_list[2], end = line_list[3],
                                                      comment = line_list[-1]))

        details_file.add_output_line(combined.get_m8()) # DirectiveDetails(line_list[0:2] + line_list[-1]))

        directives_file.write_output()
        details_file.write_output()


def main():
    obj_list = []
    if options.gcontam or (options.mito and options.gcontammitodb):
        db=constant.BLAST_MITOGCONTAM if options.gcontammitodb else constant.BLAST_GCONTAM 
        command("run_blast.py --ncbi -t " + str(options.threads) + " --output " + options.output + ".ncbi_screen.xml " +
                db + " " + options.input)
        if options.gcontam:
            gcontam = gcontam_filters(options.output + ".ncbi_screen.xml")
            obj_list.append(gcontam)
        
        if options.mito:
            mito = mito_filters(options.output + ".ncbi_screen.xml")
            obj_list.append(mito)
        
    if options.mito and not options.gcontammitodb:
        command("run_blast.py --mito -t " + str(options.threads) + " --output " + options.output + ".mito_screen.xml " +
                constant.BLAST_MITO + " " + options.input)
        ## Screen for mito
        if options.mito:
            mito = mito_filters(options.output + ".mito_screen.xml")
            obj_list.append(mito)

    ## Screen for Vector
    if options.vecscreen:
        db = constant.BLAST_ADAPTORS if options.adaptorsdb else constant.BLAST_UNIVEC
        command("run_blast.py --vecscreen -t " + str(options.threads) + " --output " + options.output + ".vecscreen.xml " +
                db + " " + options.input)
        vec = vecscreen_filters(options.output + ".vecscreen.xml")
        obj_list.append(vec)

    ## Screen for rRNA
    if options.rRNA:
        command("run_blast.py --rRNA -t " + str(options.threads) + " --output " + options.output + ".rRNAscreen.xml " +
                constant.BLAST_rRNA + " " + options.input)
        rRNA = rRNA_filters(options.output + ".rRNAscreen.xml")
        obj_list.append(rRNA)
    ## Screen for common prok contaminants
    if options.common_prok:
        db = constant.BLAST_COMMON_PROK
        cmd = "run_blast.py -b megablast -t {} --output {} {} {} ".format(str(options.threads),
                                                                             "{}.common.xml".format(options.output),
                                                                             db,
                                                                          options.input)
        command(cmd)
        common = no_filter(options.output + ".common.xml")
        obj_list.append(common)
        #command()
    ## Screen for common euk contaminants
    if options.common_euk:
        db = constant.BLAST_COMMON_EUK
        cmd = "run_blast.py -b megablast -t {} --output {} {} {} ".format(str(options.threads),
                                                                          "{}.common.xml".format(options.output),
                                                                          db,
                                                                          options.input)
        command(cmd)
        common = no_filter(options.output + ".common.xml")
        obj_list.append(common)
    ## Screen custom db
    if options.custom_db:
        command("run_blast.py -b megablast -t " + str(options.threads) + " -d --output " + options.output +
                ".custom_screen.xml " + options.custom_db + " " + options.input)
        if options.filter_as == "gcontam":
            custom = gcontam_filters(options.output + ".custom_screen.xml")
        elif options.filter_as == "mito":
            custom = mito_filters(options.output + ".custom_screen.xml")
        elif options.filter_as == "vecscreen":
            custom = vecscreen_filters(options.output + ".custom_screen.xml")
        elif options.filter_as == "rRNA":
            custom = rRNA_filters(options.output + ".custom_screen.xml")
        elif not options.filter_as:
            custom = no_filter(options.output + ".custom_screen.xml")
        else:
            print "Error filter_as option: " + options.filter_as + " is not an accepted value."
            sys.exit(-1)
        obj_list.append(custom)


    ## Combine hits, remove redundant hits, merge and extend hits
    if len(obj_list) > 1:
        combined = CombinedContaminationObject(obj_list)
    else:
        combined = obj_list[0]

    if options.redundant:
        combined.remove_redundant_hits()
        log_filter(combined, "Post Redundancy")

    if options.extend:
        combined.merge_and_extend_hits(options.reach)
        log_filter(combined, "Post Merging and Extension")

    if options.output_format.lower() == "m8":
        final_string = combined.get_m8()
    elif options.output_format.lower() == "subprep":
        final_string = list_to_string(combined.get_removal_coords(True))
    elif options.output_format.lower() == "general":
        final_string = list_to_string(combined.get_removal_coords(False))
    elif options.output_format.lower() == "directive":
        __write_directives_file(combined)
        return 0
    else:
        print "Error: Output Format can not be None"
        sys.exit(-1)

    ###Write out contamination
    if options.stdout:
        print "HIT LIST:\n" + final_string
    else:
        contam_output = open(options.output + ".contam.list", 'w')
        contam_output.write(final_string)
        contam_output.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())

