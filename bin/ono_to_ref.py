#!/usr/bin/env python

import sys
import re
import os
from gaemr.RunCommand import RunCommand
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from gaemr.PostProcessingFile import PostProcessingFile
from gaemr.DirectiveDetails import DirectiveDetails

import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()

parser = OptionParser(usage="usage: %prog [options] nucmer.delta")

parser.add_option('--output', '-o', action="store", default="tiling_ono", type='string', dest="output",
                  help='Output prefix.')
parser.add_option('--nounaligned', '-n', action="store_true", default=False, dest="nounaligned",
                  help='Do not print unaligned contigs at the end of the file.')
parser.add_option('--circular', '-c', action="store_true", default=False, dest="circular",
                  help='Assume reference sequences are circular, and cut at reference start when spanning.')
parser.add_option('--rename', '-r', action="store_true", default=False, dest="rename",
                  help='Rename ono\'d contigs and scaffolds (default=%default)')
parser.add_option('--identity', '-i', action="store", default="0.0", dest="id",
		  help='Minimum identity of alignment for inclusion in order and orienting.')
parser.add_option('--coverage', '-v', action="store", default="0.0", dest="coverage",
		  help='Minimum contig coverage to tile.')

(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("Must supply a delta file from a NUCMER alignment.")
    sys.exit(-1)

def _build_showtiling_command(delta, arg_list):
    """
    Given a delta file and argument list, return a list
    of the showtiling command and its arguments.
    """
    return [constant.SHOWTILING] + arg_list + [delta]

def _get_query_and_reference_from_delta_file(delta):
    """
    Given a nucmer delta file, read the first line to get
    the reference and query used to generate alignment.
    """
    delta_file = open(delta,'r')
    reference,query = delta_file.readline().rstrip().split(" ")
    delta_file.close()

    return reference, query

def _get_query_sequences(query):
    """
    Given a query fasta filename, open it and return a dictionary with the key
    of sequence id and the value of the sequence itself.
    """
    query_seqs = {}
    
    for scaffold in SeqIO.parse(open(query, 'r'), "fasta"):
        query_seqs[str(scaffold.id)] = str(scaffold.seq)
    
    return query_seqs

def _get_reference_details(reference):
    """
    Given a reference fasta file name, open it and return a list of the ids and 
    a dictionary of sequence lengths.
    """

    ref_list = []
    ref_lengths = {}

    for ref in SeqIO.parse(open(reference, 'r'), "fasta"):
        ref_list.append(ref.id)
        ref_lengths[ref.id] = len(ref.seq)
    
    return ref_list, ref_lengths

def _get_circular_show_tilings_command(output_header, delta):
    """
    Build up the circular show-tilings command from output header and 
    the nucmer delta file.
    """
    arg_list = ['-R', '-a', '-v 5', '-g -1', '-V 0', '-u', 
    output_header + '.unplaced','-c']

    return RunCommand(_build_showtiling_command(delta, arg_list))

def _get_circular_alignments(command):
    """
    Given a RunCommand object, return
    a dictionary of the query sequences that show circularity.
    """
    print "Executing:", command.get_command()
    print 

    circular_query_dict = {}
    for line in (command.run_command().splitlines()):
        align = line.rstrip().split("\t")
        
        # our reference alignment starts at 1, but our 
        # query alignment does not
        if int(align[0]) == 1 and not int(align[2]) == 1:
            # we only need the id and the position to cut at 
            circular_query_dict[align[12]] = int(align[2])

    return circular_query_dict

def _get_show_tilings_command(output_header, id, coverage, circular, delta):
    """
    Build up the basic show-tilings command from output header, an identity cutoff, 
    a coverage cutoff, whether or not we want circular alignments, and the nucmer
    delta file.
    """
    arg_list = ['-R', '-u', output_header + '.unplaced', "-i", id, "-v", coverage]
    if circular:
        arg_list += ["-c"]
        
    return RunCommand(_build_showtiling_command(delta, arg_list))

def _get_tilings_information(command):
    """
    Given a RunCommand object, return a dictionary
    of all the alignments (tilings).
    """
    print "Executing:", command.get_command()
    print 
    
    ref_tilings_dict = {}
    ref = ''

    for line in command.run_command().splitlines():
        if (line.startswith(">")):
            ref, bases, null = line.rstrip().split(" ")
            ref = re.sub('>','', ref)
            ref_tilings_dict[ref] = []
        else:
            ref_tilings_dict[ref].append(line.rstrip().split('\t'))
    
    return ref_tilings_dict

def _get_circular_seq_record(query_seq, query_name, cut_site, orientation, id, cov):
    """
    Given a query sequence, its name, a base position to cut at, the orientation, identity, and
    coverage of the alignment, return BioPython SeqRecords of the main sequence to print, the cut 
    sequence to print at the end of the reference, and some text describing how the cuts line up.
    """
    cut_text = ""
    if orientation == '+':
        # We already align 5'->3', so we need to get the specific subsequences (i.e. the sequence
        # that aligns starting base 1 of the reference entry, and the sequence prior to this start)
        print "\tWriting right side of " + query_name + " (" + str(cut_site) + " to end) aligned " + \
            cov + " at " + id + "% identity starting at 1 in the forward orientation."
        cut_text = "\tWriting left side of " + query_name + " (1 to " + str(cut_site) + ") aligned " + \
            cov + " at " + id + "% identity in the forward orientation."
        
        return (SeqRecord(Seq(query_seq[cut_site - 1:len(query_seq)]), query_name + '_right'),
                SeqRecord(Seq(query_seq[0:cut_site - 1]), query_name + '_left'), 
                cut_text)
    else:
        print "\tWriting right side of " + query_name + " (1 to " + str(cut_site) + ") aligned " + \
            cov, "at", id + "% identity starting at 1 in the reverse orientation."
        cut_text = "\tWriting left side of " + query_name + " (" + str(cut_site + 1) +" to end) aligned " + \
            cov + " at " + id + "% identity in the reverse orientation."
        
        return (SeqRecord(Seq(reverse_complement(query_seq[0:cut_site])), query_name + '_right'),
                SeqRecord(Seq(reverse_complement(query_seq[cut_site:len(query_seq)])), query_name + '_left'), 
                cut_text)
            
def _parse_tilings(reference_list, reference_lengths, reference_tilings, query_seqs, circular_queries):
    """
    Given a list of reference ids, their sequence lengths, the alignments or tilings, the query sequences,
    and any circular queries, return a list of ono'd BioPython SeqRecords and a list of details about the
    alignments.
    """
    ono_sequences = []
    ono_details = []

    comment = ""
    for ref in reference_list:
        print "Ordering to", ref
        cut_seq_record = None
        cut_text = ""
        
        # bail out if there is nothing to check
        if ref not in reference_tilings or not reference_tilings[ref]:
            print "Found reference entry without any alignments.  Ignoring..."
            continue

        for tiling in reference_tilings[ref]:
            rstart, rend, qstart, qend, cov, id, orientation, query = tiling

            # look for circularity first
            if int(rstart) < 0 and options.circular:
                if query in circular_queries:
                    print "\tStoring left side of the query to print at the end of alignments to this reference sequence"
                    seq_record, cut_seq_record, cut_text = _get_circular_seq_record(query_seqs[query], query, 
                                                                                    circular_queries[query], orientation, id, cov)
                    ono_sequences.append(seq_record)
                    comment = "Circular overlap found."
                    del query_seqs[query]
                    ono_details.append([ref, reference_lengths[ref], query, qend, orientation, rstart, rend, qstart, cov, id, comment])
                    continue
                else:
                    print "\tIt appears this query sequence overlaps, but no clean cut site was found.  Printing as normal."                    

            if orientation == "+":
                print "\tWriting", query, "aligned", cov, "at", id + "% identity starting at", rstart, "in the forward orientation."
                ono_sequences.append(SeqRecord(Seq(query_seqs[query]), query))
            else:
                print "\tWriting", query, "aligned", cov, "at", id + "% identity starting at", rstart, "in the reverse orientation."
                ono_sequences.append(SeqRecord(Seq(reverse_complement(query_seqs[query])), query))
    
            comment = "Alignment to reference found."
            del query_seqs[query]
            ono_details.append([ref, reference_lengths[ref], query, qend, orientation, rstart, rend, qstart, cov, id, comment])
        
        # we had a wrap that needs printing out after going through the tiling information
        if cut_seq_record:
            ono_sequences.append(cut_seq_record)
            print cut_text
    
    # if user wants the unaligned, by not specifying the nounaligned option,
    # then get the remaining query sequences and print them out
    if not options.nounaligned:
        for contig in query_seqs:
            ono_sequences.append(SeqRecord(Seq(query_seqs[contig]), contig))
            ono_details.append(["NA", "NA", contig, len(query_seqs[contig]), "NA", 
                                "NA", "NA", "NA", "NA", "NA", "No alignment found"])

    return ono_sequences, ono_details

def _write_details_to_file(ono_details_list, output_file):
    """
    Given a details list and an output file, print the contents of the list into
    a GAEMR-specific details file.
    """
    details_file = PostProcessingFile(output_file)
    details_file.add_output_line(DirectiveDetails(["REFERENCE", "REF_LENGTH", "QUERY", "QUERY_LENGTH", "ORIENTATION",
                                                   "DISTANCE_TO_NEXT_QUERY", "ALIGNMENT_COVERAGE", "AVERAGE_PCT_ID", "COMMENT"]))
    
    for i in ono_details_list:
        details_file.add_output_line(DirectiveDetails(i))

    details_file.write_output()

def main():

    delta = args[0]
    
    # get information about our inputs
    reference, query = _get_query_and_reference_from_delta_file(delta)
    query_seqs = _get_query_sequences(query)
    reference_list, reference_lengths = _get_reference_details(reference)
    
    # get alignment information
    command = _get_show_tilings_command(options.output, options.id, options.coverage, 
                                        options.circular, delta)
    reference_tilings = _get_tilings_information(command)
                               
    # see if we need circular query information
    circular_queries = {}
    if options.circular:
        command = _get_circular_show_tilings_command(options.output, delta)
        circular_queries = _get_circular_alignments(command)

    print "Ordering and orienting using", delta
    print "Reference", reference
    print "Query", query
    
    # get ono information from our gathered data
    ono_sequences_list, ono_details_list = _parse_tilings(reference_list, reference_lengths, 
                                                          reference_tilings, query_seqs, circular_queries)

    # print out our output
    _write_details_to_file(ono_details_list, options.output + ".ono.details.txt")
    interim_fasta = options.output + ".interim.fasta"
    SeqIO.write(ono_sequences_list, interim_fasta, "fasta")    
    make_assembly_command = ["make_standard_assembly_files.py", "-S", interim_fasta, "-o", options.output + ".ono"]

    if options.rename:
        make_assembly_command += ['-r']

    rc = RunCommand(make_assembly_command)
    print "Executing", rc.get_command()
    rc.run_command()

    return 0

if __name__ == "__main__":
    sys.exit(main())
