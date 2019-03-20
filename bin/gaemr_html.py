#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

#--for testing on mac, delete bfr checkin
#import sys
#sys.path.append('/Users/harindra/gaemr/trunk/gaemr')
#--


from optparse import OptionParser
import sys

from gaemr.Site import Site     #uncomment this
#from Site import Site            #delete this


#start execution here
def main():
    parser = OptionParser(usage="usage: %prog [options] assembly_name gaemr_output_dir site_name(e.g. html)")
    parser.add_option('--base_url','-b', default=None, help='base url to return to after analysis file inspection is complete')
    parser.add_option('--organism_name', default=None, help='Organism name to display in page header')
    (options, args) = parser.parse_args()

    #positional arguments only
    try:
        assembly_name = args[0]
        gaemr_run_dir = args[1]
        site_to_be_built_at = args[2]

    except:
        print "this requires 3 positional arguments\n"
        print "1: the assembly name"
        print "2: the location of the gaemr output directory"
        print "3: the gaemr subdirectory where you would like the site to be built (normally 'html')"
        print "one or more of above were missing or illegal"
        print
        sys.exit(0)

    site = Site(assembly_name,gaemr_run_dir,site_to_be_built_at,options.base_url, options.organism_name)

    site.build()

    print "site generation complete."

if __name__ == "__main__":
    main()
