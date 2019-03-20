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
import re
import operator
from Bio import SeqIO
from gaemr.RunCommand import RunCommand
import gaemr.PlatformConstant as pc
from gaemr.SimpleTable import SimpleTable
from gaemr.Taxonomy import Taxonomy

constant = pc.PlatformConstant()

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] <blast tab delimited> <query fasta file>")

parser.add_option('--nodes', '-n', action="store", default=constant.BLAST_NODES, type='string', dest="nodes_db",
                  help='NCBI nodes db file (default=%default)')
parser.add_option('--names', '-a', action="store", default=constant.BLAST_NAMES, type='string', dest="names_db",
                  help='NCBI names db file (default=%default)')
parser.add_option('--db', '-d', action="store", default=constant.BLAST_NT, type='string', dest="db",
                  help='NCBI nt db local location (default=%default)')
parser.add_option('--output', '-o', action="store", default=None, type='string',dest="output",
                   help='Output file header for contig taxonomic information (default=%default)')
parser.add_option('--heatmap_data','-m',action="store", default=None, type='string', dest="hm_data",
                  help='Output file name for heatmap-like output with contigs colored according to organism (default=%default)')
parser.add_option('--no_annotate',action="store_true",default=False,dest="no_annotate",
                  help='Whether or not to print out blast m8 with annotated sequence titles (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
                  help='Whether or not to print out html table output (default=%default)')

(options, args) = parser.parse_args()

if len(args) != 2: 
    parser.error("Must supply a blast space delimited file and query fasta file.")
    sys.exit(-1)
    
class ContigLength(object):
    
    def __init__(self, contig_id, length):
        self.contig_id = contig_id
        self.length = length

    def size(self):
        return self.length

    def name(self):
        return self.contig_id

class BlastLine(object):
    
    def __init__(self, gi, blast_line):
        self.gi = gi
        self.blast_line = blast_line
 
    def __str__(self):
        return '\t'.join([str(x) for x in self.blast_line])

    def gi_number(self):
        return self.gi

    def query(self):
        return self.blast_line[0]

    def subject(self):
        return self.blast_line[1]

    def pct_id(self):
        return float(self.blast_line[2])

    def alignment_length(self):
        return int(self.blast_line[3])

    def mismatches(self):
        return self.blast_line[4]

    def gap_opens(self):
        return self.blast_line[5]

    def subject_start(self):
        return int(self.blast_line[6])

    def subject_end(self):
        return int(self.blast_line[7])

    def query_start(self):
        return int(self.blast_line[8])
    
    def query_end(self):
        return int(self.blast_line[9])

    def evalue(self):
        return self.blast_line[10]
    
    def bit_score(self):
        return float(self.blast_line[11])

    def get_blast_line_copy(self):
        return list(self.blast_line)

# private function to get lengths of queries
def _get_lengths(file):
    """Returns length of blast query sequences."""
    lengths = []
    try:
        for s in SeqIO.parse(open(file, 'r'), "fasta"):
            lengths.append(ContigLength(s.id, len(s.seq)))
    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1
    return lengths

# private function to get blast information
def _get_blast_data(file):
    """Returns blast dictionary of lists based on query id along with a
    dictionary of blast hid id to gi number."""

    blast_data = {}
   
    try:
        with open(file) as f:
            for line in f:
                blast_line = line.rstrip('\n').split('\t')
                assert len(blast_line) == 12, 'Must supply tab-delimited BLAST output.'
                blast_line[6] = int(blast_line[6])
                blast_line[11] = float(blast_line[11])
                
                gi_hit = blast_line[1].split('|')
                gi = blast_line[1]
                if len(gi_hit) > 1:
                    gi = gi_hit[1]
                    
                if blast_line[0] not in blast_data:
                    blast_data[blast_line[0]] = []  

                blast_data[blast_line[0]].append(BlastLine(gi, blast_line))

    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1
    return blast_data

# private function to get the overlapping hits and
#   concatenate them down
def _get_gi_hits_length(blast_data):
    tmp = {}

    for i in blast_data:

        tmp[i] = {}
        prev_end = 0
        prev_start = 0
        prev_gi = None

        blast_data[i] = sorted(blast_data[i], key = lambda x: (x.subject(), x.subject_start()))
                                       
        # go through the query's hits
        for d in blast_data[i]:
            
            if not d.subject():
                continue
        
            # get some vitals about the hit
            gi = d.gi_number()
            start = d.subject_start()
            end = d.subject_end()
            length = (end - start) + 1
            
            # have we seen this gi before
            if prev_gi:

                # is it the same one?
                if gi == prev_gi:

                    # check redundancy, just in case
                    if start >= prev_start and end <= prev_end:
                        continue

                    # add the new sequence hit to what we already have
                    if start < prev_end:
                        tmp[i][gi] += (end - prev_end) + 1
                        prev_end = end
                    else: # we just add in this length if no overlap
                        tmp[i][gi] += length
                        prev_start = start
                        prev_end = end

                # first time for gi, store it
                else:
                    tmp[i][gi] = length
                    prev_start = start
                    prev_end = end
            
                prev_gi = gi
            # we don't have a previous gi
            else:
                tmp[i][gi] = length
                prev_gi = gi
                prev_start = start
                prev_end = end

    return tmp

# private function to get longest hit and gi
def _get_longest_covered_hit(data, contig_gi_hit_length):
    """Returns dict of blast hits and hit lengths"""
    blast_hit = {}
    hit_len = {}

    # for each query
    for i in data:
        # get the blast hits
        tmp = contig_gi_hit_length[i]
        
        hit_len[i] = {}

        # go through hits and make the hit length values, reverse sort sizes
        if tmp:
            for k in tmp:
                hit_len[i][k] = tmp[k]
            blast_hit[i] = sorted(tmp.iteritems(), key=operator.itemgetter(1), reverse=True)[0]

    return blast_hit,hit_len

def _get_best_hit_by_score(blast_data):
    best_score = {}
    for i in blast_data:
        tmp = sorted(blast_data[i], key = lambda x: x.bit_score(), reverse=True)[0]
        
        if not re.search('\|',tmp.subject()):
            continue
    
        gi = tmp.gi_number()
        best_score[i] = (gi, int(tmp.alignment_length()))

    return best_score

def _get_coverage_helper(id, q_length, blast_dict, tax_obj, gi_tax_dict):
    tmp = []
    
    if id in blast_dict:
        covered_len = blast_dict[id][1]
        if covered_len > q_length:
            covered_len = q_length
        covered_pct = "%.2f" % ((float(covered_len)/q_length)*100)
        tmp += [covered_len, covered_pct]        
    else:
        tmp += ["0","0.0"]

    if tax_obj.have_nodes() and tax_obj.have_names() and id in blast_dict:
        tmp.append(tax_obj.get_taxonomy_string(gi_tax_dict[blast_dict[id][0]]))
    else:
        tmp.append("domain=NO_HIT")

    return tmp

# private function to get the output string for hit covered, etc.
def _get_coverage_string(lengths, blast_hit, tax_obj, gi_tax_dict, best_score, most_common, non_genomic):
    """Returns blast coverage line."""
    covered = []
    for i in lengths: 
        q_id = i.name()
        q_length = i.size()

        tmp = [q_id,q_length]
        tmp += _get_coverage_helper(q_id, q_length, blast_hit, tax_obj, gi_tax_dict)
        tmp += _get_coverage_helper(q_id, q_length, best_score, tax_obj, gi_tax_dict)

        if q_id in most_common:
            tmp += [most_common[q_id][1]]
        else:
            tmp += ["0"]

        if q_id in non_genomic:
            tmp += [non_genomic[q_id]]
        else:
            tmp += ['None']

        covered.append(tmp)

    return covered

# public function to check for presence of dbs
def check_dbs(nodes, names, db):
    """Checks whether or not we have some dbs to use for our analysis."""
    if not nodes and not names or nodes and not names or names and not nodes:
        if not constant.BLAST_NODES or not constant.BLAST_NAMES:
            print "If using local BLAST taxonomy dump files, must have both nodes and names files as input."
            sys.exit(-1)
        nodes = constant.BLAST_NODES
        names = constant.BLAST_NAMES
        
    if not db:
        if not constant.BLAST_NT:
            print "Must supply a path to local BLAST nt database."
            sys.exit(-1)
        db = constant.BLAST_NT

    if options.hm_data and not names or not nodes:
        print "Need nodes and names file for heatmap data."
        sys.exit(-1)

    return nodes, names, db

# private function to get a heatmap-like output
def _get_heatmap_string(contig, genus_dict, lengths):
    """Returns heatmap-like string for painting a contig to taxonomic value."""
    tmp_string = "#KEY:  "
    for k,v in sorted(genus_dict.items(), key=operator.itemgetter(1)):
        tmp_string += "%s=%d;" % (k,v)
    key_string = tmp_string[:-1]
    key_string += "\n"
    
    h_string = ""
    for i in lengths:
        contig_name = i.name()
        h_string += "%s\t%s\t" % (contig_name, i.size())
        for j in xrange(int(len(contig[contig_name]))):
            h_string += "%.2f\t" % float(contig[contig_name][j])
        h_string.rstrip('\t')
        h_string += "\n"
    return key_string + h_string

# private function to get the heatmap data
def _get_heatmap_data(gi_tax_dict, tax_obj, lengths, data, hit_lengths):
    contig = {}
    genus_dict = {}
    genus_dict["No Hit"] = 0
    genus_count = 1

    for i in lengths:
        contig_name = i.name()
        contig_length = i.size()
        contig[contig_name] = [0] * 100
        tmp = {}
        for j in data[contig_name]:
            gi = j.gi_number()
            if gi not in tmp:
                tmp[gi] = []
            tmp[gi].append([int(j.subject_start()),int(j.subject_end()),int(re.sub("\.\d+","",str(j.pct_id())))])
        hits = sorted(hit_lengths[contig_name].iteritems(), key=operator.itemgetter(1), reverse=True)
        gis = []

        for j in hits:
            gis.append(j[0])

        for k in gis:
            for t_list in tmp[k]:
                if k not in gi_tax_dict:
                    continue
                genus = tax_obj.get_genus(gi_tax_dict[k])
                if genus == 'root':
                    genus = tax_obj.get_species(gi_tax_dict[k])
		if genus not in genus_dict:
                    genus_dict[genus] = genus_count
                    genus_count += 1

                start = int((float(t_list[0])/contig_length)*100)
                end = int((float(t_list[1])/contig_length)*100)
                id = "00"
                if not re.match("100",str(t_list[2])):
                    id = t_list[2]

                for num in xrange(start,end):
                    if contig[contig_name][num] == 0:
                        contig[contig_name][num] = str(genus_dict[genus]) + "." + str(id)

    return _get_heatmap_string(contig,genus_dict,lengths)

def _get_most_common_genus(data, gi_name_lookup):
    tmp = {}
    for i in data:
        for j in data[i]:
            gi = j.gi_number()
            if not gi:
                continue
            
            if gi_name_lookup[gi] not in tmp:
                tmp[gi_name_lookup[gi]] = 0
            tmp[gi_name_lookup[gi]] += int(j.alignment_length())
    
    return sorted(tmp.iteritems(), key=operator.itemgetter(1), reverse=True)[0]

def _get_most_common_genus_hits(blast_data, gi_name_lookup, genus_string, contig_gi_hit_lengths):
    common_genus_lengths = {}

    tmp_hits = {}
    for i in blast_data:
        tmp_hits[i] = []

        for j in blast_data[i]:
            if not re.search('\|', j.subject()):
                continue
            
            blast_line_copy = j.get_blast_line_copy()
            genus = gi_name_lookup[j.gi_number()]
            blast_line_copy[1] = genus
            j_copy = BlastLine(genus, blast_line_copy)
            tmp_hits[i].append(j_copy)
            #print j_copy.subject()
            
    tmp = _get_gi_hits_length(tmp_hits)

    #print tmp
    for i in tmp:
        for j in tmp[i]:
            if j == genus_string:
                common_genus_lengths[i] = (j,tmp[i][j])

    return common_genus_lengths

def _find_non_genomic_types(data, tax_obj):
    non_genomic = {}

    for i in data:

        non_genomic[i] = ' '                       
        type_dict = {}

        for j in data[i]:
            if not re.search('\|',j.subject()):
                continue

            gi = j.gi_number()
            seq_title = tax_obj.get_sequence_title_from_gi(gi)

            categorized = 0
            if re.search('(([Pp]la|[FfCc]o)smid)', seq_title):
                type_dict['PL'] = 1
                categorized = 1
            if re.search('[Vv]ector', seq_title):
                type_dict['VE'] = 1
                categorized = 1
            if re.search('[Pp]hage', seq_title):
                type_dict['PH'] = 1
                categorized = 1
            if re.search('[Vv]irus', seq_title):
                type_dict['VI'] = 1
                categorized = 1
            if re.search('[Mm]ito', seq_title) and not re.search('[Mm]itosis', seq_title):
                type_dict['MT'] = 1
                categorized = 1

            if not categorized:  #non_genomic[i] == ' ':
                type_dict['GE'] = 1

        if not type_dict:
            non_genomic[i] = 'NA'
        else:
            non_genomic[i] += ','.join([k for k in sorted(type_dict.keys())])

    return non_genomic

def _make_annotated_blast_file(blast_file, lengths, blast_data, tax_obj):

    annotated_file = re.sub("txt$","annotated.txt", blast_file)

    try:
        outfile = open(annotated_file, 'w')
    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1        

    for i in lengths:

        q_id = i.name()
        q_length = i.size()

        for d in blast_data[q_id]:
            annotation = 'NO_HIT'
            if d.subject():
                annotation = tax_obj.get_sequence_title_from_gi(d.gi_number())
            
            out_string = '\t'.join([str(d), annotation])

            outfile.write(out_string + '\n')

    outfile.close()

def _get_gi_set(blast_data):
    
    gi_list = []
    for d in blast_data:
        gi_list += [x.gi_number() for x in blast_data[d] if x.gi_number()]

    return set(gi_list)

def main():

    nodes_db, names_db, blast_db = check_dbs(options.nodes_db,options.names_db,options.db)
    
    blast_data = _get_blast_data(args[0])
    # dictionary of list of BlastLines
    
    gi_list = list(_get_gi_set(blast_data))
    
    lengths = _get_lengths(args[1])
    # lengths will be a list of ContigLength objects

    tax_obj = Taxonomy(nodes_db=nodes_db, names_db=names_db, blast_db=blast_db, gi_list=gi_list)

    
    if not options.no_annotate:
        _make_annotated_blast_file(args[0], lengths, blast_data, tax_obj)
    
    contig_gi_hit_lengths_dict = _get_gi_hits_length(blast_data)

    gi_tax_dict = tax_obj.get_gi_tax_lookup(gi_list)

    gi_name_lookup = tax_obj.get_gi_tax_name_lookup_by_level(gi_tax_dict, 'genus')
    #print contig_gi_hit_lengths_dict
    
    most_common_genus = _get_most_common_genus(blast_data, gi_name_lookup)
    #print most_common_genus

    blast_hit,hit_lengths = _get_longest_covered_hit(blast_data, contig_gi_hit_lengths_dict)

    #print blast_hit, hit_lengths

    blast_score = _get_best_hit_by_score(blast_data)
    #print blast_score

    most_common_genus_hits = _get_most_common_genus_hits(blast_data, gi_name_lookup, most_common_genus[0], contig_gi_hit_lengths_dict)
    #print most_common_genus_hits

    non_genomic_hits = _find_non_genomic_types(blast_data, tax_obj)
    
    title = "Taxonomic Classification of BLAST Hits"
    headers = ["QueryId","QueryLen","LongestHitLen","LongestPctCovered","LongestHitTaxonomy","BestScoringHitLen","BestScoringPctCovered",
               "BestScoringTaxonomy", "MostCommonOrg (" + most_common_genus[0] +")", "SequenceAnnotations"]

    t_data = _get_coverage_string(lengths, blast_hit, tax_obj, gi_tax_dict, blast_score, most_common_genus_hits, non_genomic_hits)

    st = SimpleTable(headers,t_data,title)

    output = None
    if options.output:
        output = options.output + ".blast_hit_taxonomy"
    st.print_output(output,options.html)
   
    if options.hm_data:
        try:
            output = open(options.hm_data,'w')
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            return -1            
        output.write(_get_heatmap_data(gi_tax_dict, tax_obj, lengths, blast_data, hit_lengths))
        output.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())

