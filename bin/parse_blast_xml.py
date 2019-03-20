#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import sys
from gaemr.BlastFilter import BlastFilter

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] <blast xml file>")

parser.add_option('--output', '-o', action="store", default=None, type='string', dest="output",
                  help='Parsed blast output file name (default=%default)')

parser.add_option('--print_hit_definitions', '-p', default=False,action="store_true", dest="want_hit_defs",
                  help='Print out hit definitions instead of ids in m8 output (default=%default)')

parser.add_option('--no_filter', '-n', default=False, action="store_true", dest="no_filter",
                  help='Do not filter out redundant blast hits (default=%default)')

parser.add_option('--no_hits', '-i', default=False, action="store_true", dest="no_hits",
                  help='Do not include contigs without any blast hits (default=%default)')

(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("Must supply blast xml file.")
    sys.exit(-1)


def main():

    bf = None
    if options.want_hit_defs:
        bf = BlastFilter(args[0],False)
    else:
        bf = BlastFilter(args[0])

    if not options.no_filter:
        bf.remove_redundant_hits()

    no_hits = 1
    if options.no_hits:
        no_hits = None

    bf.get_m8(no_hits)        
    m8_string = bf.get_m8(1 if not options.no_hits else None)
    
    # prepare output
    table = sys.stdout
    if options.output:
        table = open(options.output,'wb')
        sys.stdout = table

    print m8_string

    return 0

if __name__ == "__main__":
    sys.exit(main())
