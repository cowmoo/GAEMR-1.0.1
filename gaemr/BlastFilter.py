#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
import operator
import sys
import copy
from Bio.Blast import NCBIXML
from gaemr.PlatformConstant import PlatformConstant as pc



class BlastFilter(object):
    """ This class represents a filter for a BLAST alignment."""

    def __init__(self, xml_file, want_gi_ids=True, restorable=False):
        self.want_gi_ids = want_gi_ids
        self.restorable=restorable
        if xml_file:
            blast_objects = self._get_blast_obj_from_xml(xml_file)
            #     self.id_dict = self._build_ID_lookup(blast_objects)
            (self.blast_dict, self.nohit_dict, self.id_dict, self.query_dict) = self._format_blast_data(blast_objects)
            if self.restorable:
                self._orig_blast_dict=copy.deepcopy(self.blast_dict)
            
    def restore(self):
        """resore the original state, might be nice for testing purposes"""
        if self.restorable:
            self.blast_dict=copy.deepcopy(self._orig_blast_dict)

    # private function to parse xml to Blast object
    def _get_blast_obj_from_xml(self, xml_file):
        blast_records = NCBIXML.parse(open(xml_file))
        return blast_records

    # private function to calculate % id
    def _calculate_identity(self, identities, align_length):
        return round(((float(identities) / align_length) * 100), 2)

    # private function to make a blast m8-like line from the hit
    def _make_blast_list(self, query, hit_id, hsp):
        """
        returns 12 element list 
        0: query (string)
        1: hit_id (string)
        2: pct_identity (float)
        3: align_length (int)
        4: mismatches (int)
        5: gaps (int)
        6: query_start (int)
        7: query_end (int)
        8: sbjct_start (int) 
        9: sbjct_end (int)
        10: expect (float)
        11: score (flost)
        """ 
        
        # Length of the alignment
        # E.g. 2461
        # Xml: <Hsp_align-len>2461</Hsp_align-len>
        align_length=hsp.align_length
        
        # Number of identities
        # E.g. 2439
        # Xml: <Hsp_identity>2439</Hsp_identity>
        identities=hsp.identities
        
        # Number of gaps
        # E.g.: 14
        # Xml: <Hsp_gaps>14</Hsp_gaps>
        gaps=hsp.gaps
        
        # The start residue for the query sequence.  (1-based)
        # E.g.: 1
        # Xml: <Hsp_query-from>1</Hsp_query-from>
        query_start=hsp.query_start
        
        # The end residue for the query sequence.  (1-based)
        # E.g.: 2456
        # Xml: <Hsp_query-to>2456</Hsp_query-to>
        query_end=hsp.query_end
        
        # The start residue for the sbjct sequence.  (1-based)
        # E.g.: 287428
        # Xml: <Hsp_hit-from>287428</Hsp_hit-from>
        sbjct_start=hsp.sbjct_start

        # The end residue for the sbjct sequence.  (1-based)
        # E.g.: 284977
        # Xml: <Hsp_hit-to>284977</Hsp_hit-to>
        sbjct_end=hsp.sbjct_end

        # Expect value.  (float)
        # E.g.: 1.2675219251051e-09
        # Xml: <Hsp_evalue>1.2675219251051e-09</Hsp_evalue>
        expect=hsp.expect
        
        # BLAST score of hit.  (float)
        # E.g.: 2333
        # Xml: <Hsp_score>2333</Hsp_score>
        score=hsp.score
        
        # Note: Corrected this calculation because Blast does not count gaps as mismatches
        # mismatches = align_length - hsp.identities
        mismatches = align_length - identities - gaps

        pct_identity = self._calculate_identity(identities, align_length)
        return [query, hit_id, pct_identity, align_length, mismatches, gaps, query_start,
                query_end, sbjct_start, sbjct_end, expect, score]

    # private function to get the blast hits from xml
    def _get_hits(self, query, hit_id, hsp, hits):
        hits[query]['hits'].append(self._make_blast_list(query, hit_id, hsp))
        return hits

    # private function to format the blast object into our data structure
    def _format_blast_data(self, blast_records):
        hits = {}
        nohit = {}
        hit_id_lookup = {}
        query_id_lookup = {}
        for blast_record in blast_records:
            # Name of query sequence
            # E.g.: contig000761
            # Xml: <Iteration_query-def>contig000761 this contig rocks</Iteration_query-def>
            query = blast_record.query.split(' ')[0]
            
            # Length of the query
            # E.g.: 2603
            # Xml: <Iteration_query-len>2603</Iteration_query-len>
            query_length = blast_record.query_length
            
            # Identifier of the query sequence
            # Query_1
            # Xml: <Iteration_query-ID>Query_1</Iteration_query-ID>
            query_id = blast_record.query_id
            
            query_id_lookup[query] = query_id
            hits[query] = {}
            hits[query]['length'] = query_length
            hits[query]['hits'] = []
            
            # A list of Alignment objects
            # Xml: <Iteration_hits/>
            # each Alignment object holds information about one alignment hit
            # Xml: <Hit/>
            
            blast_alignments = blast_record.alignments
            
            if len(blast_alignments) == 0:
                nohit[query] = 1
                hits[query]['hits'].append([query, "", float(0), 0, 0, 0, 0, 0, 0, 0, float(0), float(0)])
                continue

            for alignment in blast_alignments:
            #   id_lookup[self._strip_pipe(alignment.hit_id)] = alignment.hit_def
                # Hit identifier
                # E.g.: gi|331223848|ref|NW_003526561.1|
                # Xml: <Hit_id>gi|331223848|ref|NW_003526561.1|</Hit_id>
                hit_id=alignment.hit_id
                
                # Hit definition
                hit_def=alignment.hit_def
                # E.g.: <Hit_def>Puccinia graminis
                # Xml: <Hit_def>w<Hit_def>Puccinia graminis</Hit_def>
                
                hit_id_lookup[hit_id] = hit_def
                
                # A list of HSP objects
                # Xml: <Hit_hsps/>
                # A HSP stores information about one hsp in an alignment hit
                # Xml: <Hsp/>
                hsps = alignment.hsps
                for hsp in hsps:
                    if self.want_gi_ids:
                        hits = self._get_hits(query, hit_id, hsp, hits)
                    else:
                        hits = self._get_hits(query, hit_def, hsp, hits)

        return hits, nohit, hit_id_lookup, query_id_lookup

    def update_nohit_list(self):
        obj = self.blast_dict
        list = self.nohit_dict
        for query in obj:
            if not query in list:
                if len(obj[query]['hits']) == 0:
                    list[query] = 1

    # public function to remove queries without hits from object
    def remove_queries_with_no_hits(self):
        obj = self.blast_dict
        list = self.nohit_dict
        copy = obj.copy()
        for query in copy:
            if query in list or len(copy[query]['hits']) == 0:
                del obj[query]
        # this should probably be cleared at this point so that 
        # the methods can still use this dictionary to get an accurate value
        self.nohit_dict.clear()

    # public function to remove completely redundant hits from object
    def remove_redundant_hits(self):
        obj = self.blast_dict
        list = self.nohit_dict
        for query in obj:
            if not query in list:
                #print obj[query]['hits']
                tmp_sort = sorted(obj[query]['hits'], key=operator.itemgetter(0, 7), reverse=True)
                final_sort = sorted(tmp_sort, key=operator.itemgetter(0, 6))

                start = 0
                end = 0
                obj[query]['hits'] = final_sort
                copy = obj[query]['hits'][:]
                for hit in copy:
                    if start <= hit[6] and end >= hit[7]:
                        obj[query]['hits'].remove(hit)
                        continue
                    start = hit[6]
                    end = hit[7]

    #private function to extend hits to contig edge if <= min_contig_size
    def _extend_hit(self, hit, length, min_size=200):
        if hit[6] <= min_size:
            #print "Extending " + hit[0] + " : " + str(hit[6]) + " to 0"
            hit[6] = 0

        if length - hit[7] <= min_size:
            #print "Extending " + hit[0] + " : " + str(hit[7]) + " to " + str(length)
            hit[7] = length

        return hit


    #public function to merge overlapping hits
    def merge_and_extend_hits(self, min_size=1):
        obj = self.ignore_hits_by_dict(self.nohit_dict)
        extended = 0
        for query in obj:
            merged_list = []
            length = obj[query]['length']
            if len(obj[query]['hits']) > 1: #has more than one hit
                tmp_sort = sorted(obj[query]['hits'], key=operator.itemgetter(0, 7), reverse=True)
                final_sort = sorted(tmp_sort, key=operator.itemgetter(0, 6))
                obj[query]['hits'] = final_sort
                merged_hit = self._extend_hit(obj[query]['hits'][0], length, min_size)

                #Create unique ID to store the merged names of the hits
                merged_id = str(merged_hit[0]) + "_" + str(merged_hit[6])
                self.id_dict[merged_id] = self._get_description_for_hit(merged_hit) # self.id_dict[merged_hit[1]]
                merged_hit[1] = merged_id

                copy = obj[query]['hits'][1:]
                for hit in copy:
                    hit = self._extend_hit(hit, length, min_size)

                    if merged_hit[7] + min_size >= hit[6]: #hits overlap so reset end coordinate
                        merged_hit[7] = hit[7]

                        #If hit name is not in the merged list add it
                        if not str(self.id_dict[merged_id]).count(str(self.id_dict[hit[1]])):
                            self.id_dict[merged_id] += ", " + self._get_description_for_hit(hit) # self.id_dict[hit[1]]
                        extended += 1

                    else: # hits don't overlap so move on to next hit
                        merged_list.append(merged_hit)
                        merged_hit = hit
                        merged_id = str(merged_hit[0]) + "_" + str(merged_hit[6])
                        self.id_dict[merged_id] = self._get_description_for_hit(merged_hit) # self.id_dict[merged_hit[1]]
                        merged_hit[1] = merged_id

                merged_list.append(merged_hit)
                obj[query]['hits'] = merged_list

            elif len(obj[query]['hits']): # Only has one hit just extend to contig ends but don't try to merge
                obj[query]['hits'][0] = self._extend_hit(obj[query]['hits'][0], length, min_size)
        print "Merged and Extended " + str(extended) + " Hits"


    #public function to filter queries from the object based on a dict
    def filter_hits_by_dict(self, in_dict, keep):
        obj = self.blast_dict
        filtered_hits = {}
        for query in obj:
            if keep:
                if query in in_dict:
                    filtered_hits[query] = obj[query]
            else:
                if not query in in_dict:
                    filtered_hits[query] = obj[query]
        return filtered_hits

    #public function to remove queries from the object based on a dict
    def ignore_hits_by_dict(self, in_dict):
        return self.filter_hits_by_dict(in_dict, False)

    #public function to select queries from the object based on a dict
    def select_hits_by_dict(self, in_dict):
        return self.filter_hits_by_dict(in_dict, True)

    # private function to determine what to return
    def _decider(self, choice=None):
        if choice == 1:
            return self.blast_dict
        else:
            return self.ignore_hits_by_dict(self.nohit_dict)

    # public function to print full blast object
    def print_object(self, unaligned=None):
        obj = self._decider(unaligned)
        for query in obj:
            print query
            print "\tLength: " + str(obj[query]['length'])
            print "\tHits:"
            for hit in obj[query]['hits']:
                #  print "\t" + self.id_dict[self._strip_pipe(hit[1])]
                if hit[1] in self.id_dict:
                    print "\tID:\t" + self.id_dict[hit[1]]
                print "\tAlign:" + ''.join([str(i) + "\t" for i in hit])

    #public function to print hits in the m8 format
    def get_m8(self, unaligned=None):
        m8_string = ""
        hit_strings = []
        hits = self._decider(unaligned)
        for query in hits:
            for hit in hits[query]['hits']:
                if isinstance(self, VecScreenFilter) or isinstance(self, MitoFilter) or isinstance(self, gContamFilter) or isinstance(self, CombinedContaminationObject):
                    hit_desc = self._get_description_for_hit(hit)
                else:
                    hit_desc = hit[1]
                # fixed here so this method does not alter state of object
                items=hit[0:1]+[hit_desc]+hit[2:]
                hit_string = "\t".join([str(i) for i in items])
                hit_strings.append(hit_string)
        m8_string = "\n".join(hit_strings)
        return m8_string

    #public function to filter by length and pct id
    def filter_by_length_and_pct_id(self, length, pct_id):
        hits = self._decider('False')
        filtered_hits = []
        for query in hits:
            #print "HIT:", hits[query]['hits']
            for hit in hits[query]['hits']:
                if int(hit[3]) >= int(length) and float(hit[2]) >= float(pct_id):
                    filtered_hits.append(hit)
        return filtered_hits

    #public function to filter by length
    def filter_by_length(self, length):
        return self.filter_by_length_and_pct_id(length, '0.0')

    #public function to filter by pct id
    def filter_by_pct_id(self, pct_id):
        return self.filter_by_length_and_pct_id(0, pct_id)

    #public function to filter hits by pct of contig covered
    def filter_by_pct_covered(self, pct_cvd):
        hits = self._decider('False')
        filtered_hits = []
        for query in hits:
            bases_covered = 0
            start = 0
            end = 0
            query_length = hits[query]['length']

            tmp_sort = sorted(hits[query]['hits'], key=operator.itemgetter(0, 7), reverse=True)
            final_sort = sorted(tmp_sort, key=operator.itemgetter(0, 6))

            for hit in final_sort:
                if start <= hit[6] and end >= hit[7]:
                    continue
                else:
                    #print hit[3]
                    bases_covered += hit[3]
                start = hit[6]
                end = hit[7]

            #print query + " - " + str(bases_covered)

            if bases_covered / query_length <= pct_cvd:
                for hit in hits[query]['hits']:
                    filtered_hits.append(hit)
        return filtered_hits

    # A bunch of new functions to help interrogate this complex object
    # might help with some unit tests
    
    #public finction to return the nubmer of hits for a particular query    
    def count_hits_for_query(self, query):
        obj=self._decider(0) 
        return len(obj[query]['hits']) if query in obj else 0

    #public function to print list of queries with no hits
    def print_queries_without_hits(self):
        for query in self.nohit_dict:
            print query
        
    #public function to return number of queries with no hits
    def count_queries_without_hits(self):
        return len(self.nohit_dict)
    
    #public function to return list of the queries with no hits
    def list_queries_without_hits(self):
        return self.nohit_dict.keys()
    
    #public function to print list of all queries with or without hits
    def print_queries_with_hits(self):
        return self._print_queries()
        
    #public function to return number of queries with with or without hits
    def count_queries_with_hits(self):
        return self._count_queries()
    
    #public function to return list of the queries with or without hits
    def list_queries_with_hits(self):
        return self._list_queries()

    #public function to print list of all queries with or without hits
    def print_queries(self):
        return self._print_queries(True)
        
    #public function to return number of queries with with or without hits
    def count_queries(self):
        return self._count_queries(True)
    
    #public function to return list of the queries with or without hits
    def list_queries(self):
        return self._list_queries(True)
    
    def _print_queries(self, without_hits=False):
        for query in self._decider(1 if without_hits else 0):
            print query

    def _count_queries(self, without_hits=0):
        return len(self._decider(1 if without_hits else 0))
    
    def _list_queries(self, without_hits=0):
        return self._decider(1 if without_hits else 0).keys()
    
    # END: A bunch of new functions to help interrogate this complex object

    #public function to print hits in the SubmissionPrep format
    def get_SubmissionPrep(self):
        self.get_removal_coords(True)

    #private function returns the 0-based id for a query
    def _get_id(self, query):
        return int(self.query_dict[query].replace("Query_", "")) - 1

    #public function to print hits in the general removal format
    ## contig_name\tcontig_len\thit_start\thit_end
    def get_removal_coords(self, as_id=False):
        out_list = []
        hits = self._decider(False)
        for query in hits:
            length = hits[query]['length']
            for hit in hits[query]['hits']:
                name = hit[0]
                if as_id:
                    name = self._get_id(hit[0])
                out_list.append([name, length, hit[6], hit[7], 0, self._get_description_for_hit(hit)])
                ### print hit
        return out_list

    # Sometimes we may use a blast database that does not have definition lines for all entries, but
    # we want a descriptive string, so the id field is the best fallback
    def _get_description_for_hit(self, hit):
        """
        If the bast hit definition is 'No definition line' (or other configurable value) returns
        the hit id instead, otherwise returns the hit definition
        """
        desc=None
        hit_id=hit[1]
        hit_def=self.id_dict[hit_id]
        
        if not hit_def or hit_def == pc.BLAST_XML_NO_DEF_LINE:
            desc=hit_id
        else:
            desc=hit_def
        
        print "DESCRIPTION: ", desc
        return desc

class rRNAFilter(BlastFilter):
    def __init__(self, xml):
        super(rRNAFilter, self).__init__(xml)


class VecScreenFilter(BlastFilter):
    """Bins hits into four catagories, String, Moderate, Weak, or No Hit based on the raw hsp score
       Different sets of minimal score parameters are applied depending on whether the hit is 
       near the end of the query contig
       """
       
    
    def __init__(self, xml, use_protection=True, protected_terms=[], proximity=25, strong=(30, 24), moderate=(25, 19), weak=(23, 16)):
        """
        use_protection: enables a special protection feature in which WEAK hits are protected from removal if the hit's id or definition 
                       contains a protected term found
        protected_terms: a list of protected terms defaults to ["illumina", "adaptor", "adapter", "titanium linker"]
        proximity: defines how close to the end a hit needs to be on the query contig to be considered a terminal hit
        strong: the score parameters for STRONG hits
        moderate: the score parameters for MODERATE hits
        weak:  the score parameters for WEAK hits
        """
        super(VecScreenFilter, self).__init__(xml)
        # protection is on by default for backwards compatability
        # old protected terms preserved as default for backwards compatability
        self.use_protection=use_protection
        self.protected_terms=protected_terms if protected_terms else ["illumina", "adaptor", "adapter", "titanium linker"]
        self.proximity=proximity
        self.strong=strong
        self.moderate=moderate
        self.weak=weak
        
    def _terminal(self, hit, query_length):
        end_range = query_length - self.proximity
        return hit[7] <= self.proximity or hit[6] <= self.proximity or hit[7] > end_range or hit[6] > end_range
    
    # private class to determine the strength of the hit
    def _get_hit_strength(self, hit, query_length):
        # index into tuple is 1 if hit is near terminus of query contig, otherwise 0
        x=1 if self._terminal(hit, query_length) else 0
        
        if self._get_score(hit) >= self.strong[x]:
            return "Strong"
        elif self._get_score(hit) >= self.moderate[x]:
            return "Moderate"
        elif self._get_score(hit) >= self.weak[x]:
            return "Weak"
        else:
            return "No Hit"


    # OLD WAY OF GETTING SCORE
    # private class to calc the score of the hit
    # hit[2] percent identity
    # hit[3] align_length
    # hit[4] mismatches
    # hit[5] gaps
    # NOTE: This calculation is not strictly the BLAST score calculation
    # BLAST scores gap openings and gap extensions separately.
    # This calculation provides an accurate value only if the alignment being scored
    # has only single base length gaps
#    def _get_score(self, hit):
#        numid = round(hit[3] * hit[2] / 100, 0)
#        return 4 * numid - 3 * hit[3] - 3 * hit[5] - 2 * hit[4]
    
    # NEW WAY OF GETTING SCORE
    # just use the raw score from the BLAST output
    def _get_score(self, hit):
        return hit[11]

    # private class to check if hit should be protected from filtering
    # NCBI determines the list of terms to be protected
    def _is_protected(self, hit):
        is_protected=False
        if self.use_protection:
            for term in self.protected_terms:
                if re.search(term, self.id_dict[hit[1]], re.IGNORECASE):
                    is_protected=True
                    break
        return is_protected

    # Remove the hits with strength of No Hit
    # Remove WEAK hits unless they are protected by the speical protection feature
    def filter_weak_hits(self):
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            for hit in copy:
                strength = self._get_hit_strength(hit, hits[query]['length'])
                if (strength == "Weak" and not self._is_protected(hit)) or strength == "No Hit" :
                    hits[query]['hits'].remove(hit)
                    #print "Ignore hit with score of ", hit[11], " and category ", strength
                    #print hit
                    continue
        self.remove_queries_with_no_hits()

    # Select the hits that contain "Illumina" in the description
    # Note: if protection is used the default protected terms should include Illumina
    def select_Illumina_hits(self):
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            for hit in copy:
                if not re.search("Illumina", self.id_dict[hit[1]]):
                    hits[query]['hits'].remove(hit)
                    continue
        self.remove_queries_with_no_hits()


class MitoFilter(BlastFilter):
    def __init__(self, xml, mito_param):
        super(MitoFilter, self).__init__(xml)
        self.mito_param=mito_param

    def select_valid_mito_hits(self):
        """old style filter which does two things:
        1) select for those hits that contain the term mitochondrion because we may use db that is a combintation
           a gcontam1 and mito.
        2) exclude hits if the aligned length of the query contig is not more than a certain percentage of the total query length
        """
        #mito_dict = {}
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            query_length = hits[query]['length']
            for hit in copy:
                if not re.search("mitochondrion", self.id_dict[hit[1]])  or  float(
                    (hit[7] - hit[6]) / float(query_length)) < self.mito_param: #or not re.search("mitochondrion",hit[1])
                    hits[query]['hits'].remove(hit)
                    continue
                covered = float((hit[7] - hit[6]) / float(query_length))
                #print self.id_dict[hit[1]] + "\t" + str(covered)
            sys.stdout.flush()
        self.remove_queries_with_no_hits()
        
    def select_hits_by_length(self):
        """remove hits if the length of the aligned query is less than cutoff value""" 
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            for hit in copy:
                if (hit[7] - hit[6] + 1 < self.mito_param):
                    hits[query]['hits'].remove(hit)
                    continue
        self.remove_queries_with_no_hits()



class gContamFilter(BlastFilter):
    
    GCONTAM_DEFAULT_CRITERIA=((98, 50),
                              (94, 100),
                              (90, 200))
    
    def __init__(self, xml, criteria=GCONTAM_DEFAULT_CRITERIA):
        super(gContamFilter, self).__init__(xml)
        self.criteria=criteria

    def ignore_mito_hits(self):
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            for hit in copy:
                if re.search("mito", self.id_dict[hit[1]]) or re.search("mito", hit[1]):
                    hits[query]['hits'].remove(hit)
                    continue
        self.remove_queries_with_no_hits()

    ## Filter hits based on the NCBI requirements for trimming
    def filter_hits(self):
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            for hit in copy:
                length = hit[3]
                pct_id = hit[2]
                if not self._match_criteria(pct_id, length):
                    hits[query]['hits'].remove(hit)
                   # print "Ignore gcontam"
                   # print hit
                   # print str(pct_id) + "\t" + str(length)
                    continue
        self.remove_queries_with_no_hits()
        
    def _match_criteria(self, pct_id, length):
        match=False
        for c in self.criteria:
            if pct_id >= c[0] and length >= c[1]:
                match=True
                break
        return match

class CombinedContaminationObject(BlastFilter):
    def __init__(self, obj_list):
        super(CombinedContaminationObject, self).__init__(None)
        self.blast_dict = {}
        self.nohit_dict = {}
        self.id_dict = {}
        self.query_dict = {}

        for obj in obj_list:
            self.id_dict.update(obj.id_dict)
            self.query_dict.update(obj.query_dict)
            if not isinstance(obj, MitoFilter) and not isinstance(obj, gContamFilter) and not isinstance(obj,
                VecScreenFilter) and not isinstance(obj, rRNAFilter):
                print "ERROR: Trying to combine non BlastFilter objects"
                sys.exit(-1)
            for query in obj.blast_dict:
                if not query in self.blast_dict:
                    self.blast_dict[query] = {}
                    self.blast_dict[query]['length'] = obj.blast_dict[query]['length']
                    self.blast_dict[query]['hits'] = []
                    self.blast_dict[query]['hits'] = obj.blast_dict[query]['hits'][:]
                else:
                    if obj.blast_dict[query]['length'] == self.blast_dict[query]['length']:
                    # print "APPEND:" + query
                    # print "ORIGINAL: ",
                    # print self.blast_dict[query]['hits']
                        self.blast_dict[query]['hits'] += obj.blast_dict[query]['hits'][:]
                        # print "APPENDED: ",
                        # print self.blast_dict[query]['hits']

                    else:
                        print "ERROR:  Trying to combine Blast filter objects from different queries: " + query + "has differing lengths " + str(
                            obj[query]['length']) + " != " + str(self.blast_dict[query]['length'])
                        sys.exit(-1)
