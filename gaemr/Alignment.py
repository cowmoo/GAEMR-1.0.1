#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import os
import subprocess
import sys
import re
from BamFile import BamFile
import PlatformConstant as pc
import logging
from os.path import join

constant = pc.PlatformConstant()


def command(cmd):
    print cmd
    subprocess.check_call(cmd, shell=True)


class Alignment(object):

    def __init__(
            self,
            query,
            reference,
            paired=True,
            short=True,
            long=False,
            mem=True,
            threads=2,
            output_header='aligned',
            output_path='.',
            force_unpaired=False,
    ):
        self.query = query
        self.reference = reference
        self.paired = paired
        self.short = short
        self.long = long
        self.mem = mem
        self.threads = threads
        self.output_path = output_path

        #    self.output_path = os.path.dirname(os.path.abspath(output_path))

        print 'header:' + output_header
        self.output_header = os.path.basename(output_header)
        self.force_unpaired = force_unpaired
        self.logger = logging.getLogger('Alignment')
        ch = logging.StreamHandler()
        fh = logging.FileHandler('Alignment.log')
        formatter = \
            logging.Formatter('[%(asctime)s | %(name)s | %(levelname)s] %(message)s'
                              )
        ch.setFormatter(formatter)
        ch.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.addHandler(fh)

    def is_paired(self):
        return self.paired

    def is_short(self):
        return self.short

    def is_long(self):
        return self.long

    def uses_mem(self):
        return self.mem

    def is_forced_unpaired(self):
        return self.force_unpaired

    def build_SeqDict(self):
        (ref, ext) = os.path.splitext(os.path.abspath(self.reference))
        if os.path.exists(ref + '.dict'):
            os.unlink(ref + '.dict')
        command('java -jar ' + constant.PICARD
                + '/picard.jar CreateSequenceDictionary  R= '
                + self.reference + ' O= ' + ref + '.dict TMP_DIR= '
                + constant.PICARD_TMP_DIR)

    def index_exists(self):
        if "/" in self.reference:
            ref_dir = "/".join(self.reference.split('/')[0:-1])
            ref_head = (self.reference.split('/')[-1]).split('.')[0]
        else:
            ref_dir = os.getcwd()
            ref_head = self.reference.split('.')[0]
        index_file_count = 0
        index_files = ['.dict', '.bwt']
        dir_contents = os.listdir(ref_dir)
        for item in dir_contents:
            for ext in index_files:
                if ref_head in item and item.endswith(ext):
                    index_file_count += 1
        if index_file_count == 2:
            return True
        else:
            return False

    def mark_duplicates(self, bam):
        output_header = bam.rstrip('.bam')
        command('java -jar ' + constant.PICARD
                + '/picard.jar MarkDuplicates I= ' + bam + ' O= '
                + output_header + '.duplicates_marked.bam M= '
                + output_header
                + '.duplicates_marked.metrics  TMP_DIR= '
                + constant.PICARD_TMP_DIR)

    def create_bam_index(self, bam):
        command(constant.SAMTOOLS + ' index ' + bam)

    def merge_fastq_files(self, list):

        output = re.sub('read1', 'merged', list[0])
        command('cat ' + list[0] + ' ' + list[1] + ' > ' + output)
        return [output]

    def convert_to_fastq(self, tmp_dir):
        if self.is_paired() or self.is_forced_unpaired():
            direction = 'fr'
        else:
            direction = None

        bam = BamFile(self.query, 'bam', tmp_dir + '/'
                      + self.output_header + '.bam', direction)
        (code, fastq_list) = bam.convert_to_fastq()
        if not code:
            print 'Conversion to Fastq failed'
        if self.is_forced_unpaired():
            fastq_list = self.merge_fastq_files(fastq_list)

        return fastq_list


class BwaAlignment(Alignment):

    def __init__(
            self,
            query,
            reference,
            output_header,
            paired=True,
            short=True,
            long=False,
            use_mem=True,
            threads=2,
            forced=False,
    ):
        full_path = os.path.dirname(os.path.abspath(output_header))
        self.logger = logging.getLogger('Alignment')
        full_header = os.path.basename(output_header)
        self.logger.info('\n=========ALIGNMENT MODULE LOGGING BEGIN========='
                         )
        self.use_mem = use_mem
        self.uses_bwasw = False
        super(BwaAlignment, self).__init__(
            query,
            reference,
            paired,
            short,
            long,
            use_mem,
            threads,
            full_header,
            full_path,
            forced,
        )

    def build_index(self, bwa_mem=False):

        if not os.path.exists(self.reference):
            print 'ERROR: Reference file: ' + self.reference \
                  + ' - Does not exist!!'
            sys.exit(-1)

        if not os.path.dirname(os.path.abspath(self.reference)) \
                == os.path.abspath(self.tmp_dir):
            ref_link = os.path.abspath(self.tmp_dir) + '/' \
                       + os.path.basename(self.reference)
            print 'Reference path does not match temp_dir - link reference to tmp dir: ' \
                  + os.path.dirname(os.path.abspath(self.reference)) \
                  + ' -TMP: ' + os.path.abspath(self.tmp_dir)
            if os.path.exists(ref_link):
                os.unlink(ref_link)
            os.symlink(os.path.abspath(self.reference), ref_link)
            self.reference = ref_link

        if os.path.getsize(self.reference) < 2000000000:
            algorithm = 'is'
        else:
            algorithm = 'bwtsw'
        if not bwa_mem:
            command(constant.BWA + '/bwa index -a ' + algorithm + ' '
                    + self.reference)
        else:
            command(constant.BWA_MEM + '/bwa index ' + self.reference)

    def __run_BWA_sam(self):
        bam = self.query
        if self.is_paired():
            command(constant.BWA + '/bwa sampe -a 100000' + ' -f '
                    + self.tmp_dir + '/' + self.output_header
                    + '_bwa.sam ' + self.reference + ' ' + self.tmp_dir
                    + '/' + self.output_header + '.1.sai '
                    + self.tmp_dir + '/' + self.output_header
                    + '.2.sai ' + bam + ' ' + bam)
        else:
            command(constant.BWA + '/bwa samse' + ' -f ' + self.tmp_dir
                    + '/' + self.output_header + '_bwa.sam '
                    + self.reference + ' ' + self.tmp_dir + '/'
                    + self.output_header + '.0.sai ' + bam)

    def __run_BWA_aln(self, mate):

        # remove -q 5 from defaults. --bruce 9/14/12

        command(constant.BWA + '/bwa aln ' + self.reference
                + ' -l 32 -k 2 -t ' + str(self.threads) + ' -o 1 -f '
                + self.tmp_dir + '/' + self.output_header + '.'
                + str(mate) + '.sai -b' + str(mate) + ' ' + self.query)

    def __run_BWA_mem(self, paired=True):
        query_fastq = self.convert_to_fastq(self.tmp_dir)

        if paired:
            command(constant.BWA_MEM + '/bwa mem ' + self.reference
                    + ' ' + query_fastq[0] + ' ' + query_fastq[1]
                    + ' > ' + self.tmp_dir + '/' + self.output_header
                    + '_bwa.sam')
        else:
            command(constant.BWA_MEM + '/bwa mem ' + self.reference
                    + ' ' + query_fastq[0] + ' > ' + self.tmp_dir + '/'
                    + self.output_header + '_bwa.sam')

    def __add_unmapped_reads(self, no_index=False, unsorted=False):
        ref_dict = self.tmp_dir + "/" + (self.reference.split('/')[-1]).split('.')[0] + ".dict"
        if not os.path.exists(ref_dict):
            print "CREATING SEQUENCE DICTIONARY..."
            self.build_SeqDict()
        paired = 'true'
        if not self.is_paired():
            paired = 'false'
        long_read_args = ''
        if not self.is_short():
            long_read_args = ' MAX_GAPS= -1 '

        output_header = ''
        if not unsorted:
            output_header = self.output_header
        else:
            output_header = self.output_header + '_unsorted'

        if self.uses_mem() and not self.uses_bwasw:
            if constant.BWA_MEM_VERSION:
                command('java -jar ' + constant.PICARD
                        + '/MergeBamAlignment.jar R= ' + self.reference
                        + ' UNMAPPED_BAM= ' + self.query + ' ALIGNED_BAM= '
                        + self.tmp_dir + '/' + self.output_header
                        + '_bwa.sam OUTPUT= ' + self.tmp_dir + '/'
                        + output_header
                        + '.bam ORIENTATIONS=FR  PROGRAM_RECORD_ID=BWA PROGRAM_GROUP_VERSION='
                        + constant.BWA_MEM_VERSION
                        + ' PROGRAM_GROUP_COMMAND_LINE=tk VALIDATION_STRINGENCY=SILENT CLIP_ADAPTERS=false '
                          'ALIGNER_PROPER_PAIR_FLAGS=true PE= '
                        + paired + ' TMP_DIR= ' + self.tmp_dir
                        + long_read_args)
            else:
                command('java -jar ' + constant.PICARD
                        + '/MergeBamAlignment.jar R= ' + self.reference
                        + ' UNMAPPED_BAM= ' + self.query + ' ALIGNED_BAM= '
                        + self.tmp_dir + '/' + self.output_header
                        + '_bwa.sam OUTPUT= ' + self.tmp_dir + '/'
                        + output_header
                        + '.bam ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CLIP_ADAPTERS=false '
                          'ALIGNER_PROPER_PAIR_FLAGS=true PE= '
                        + paired + ' TMP_DIR= ' + self.tmp_dir
                        + long_read_args)
        else:
            print "MEEP"
            if constant.BWA_VERSION:
                command('java -jar ' + constant.PICARD
                        + '/MergeBamAlignment.jar R= ' + self.reference
                        + ' UNMAPPED_BAM= ' + self.query + ' ALIGNED_BAM= '
                        + self.tmp_dir + '/' + self.output_header
                        + '_bwa.sam OUTPUT= ' + self.tmp_dir + '/'
                        + output_header
                        + '.bam ORIENTATIONS=FR  PROGRAM_RECORD_ID=BWA PROGRAM_GROUP_VERSION='
                        + constant.BWA_VERSION
                        + ' PROGRAM_GROUP_COMMAND_LINE=tk VALIDATION_STRINGENCY=SILENT CLIP_ADAPTERS=false '
                          'ALIGNER_PROPER_PAIR_FLAGS=true PE= '
                        + paired + ' TMP_DIR= ' + self.tmp_dir
                        + long_read_args)
            else:
                command('java -jar ' + constant.PICARD
                        + '/MergeBamAlignment.jar R= ' + self.reference
                        + ' UNMAPPED_BAM= ' + self.query + ' ALIGNED_BAM= '
                        + self.tmp_dir + '/' + self.output_header
                        + '_bwa.sam OUTPUT= ' + self.tmp_dir + '/'
                        + output_header
                        + '.bam ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CLIP_ADAPTERS=false '
                          'ALIGNER_PROPER_PAIR_FLAGS=true PE= '
                        + paired + ' TMP_DIR= ' + self.tmp_dir
                        + long_read_args)
        # Create BAM Index
        # self.create_bam_index(self.output_path  + "/" + self.output_header + ".bam")

    def __run_BWA_bwasw(self):
        self.uses_bwasw = True
        query_fastq = self.convert_to_fastq(self.tmp_dir)
        command(constant.BWA + '/bwa bwasw -t ' + str(self.threads)
                + ' -f ' + self.tmp_dir + '/' + self.output_header
                + '_bwa.sam ' + self.reference + ' ' + query_fastq[0])

        # Convert to BAM

        aligned_sam = BamFile(self.tmp_dir + '/' + self.output_header
                              + '_bwa.sam', 'sam', self.tmp_dir + '/'
                              + self.output_header + '_unsorted.bam')
        aligned_sam.convert_to_bam()

        # Sort BAM

        aligned_bam = BamFile(self.tmp_dir + '/' + self.output_header
                              + '_unsorted.bam', 'bam',
                              self.output_path + '/'
                              + self.output_header)
        aligned_bam.sort_bam()

        self.create_bam_index(self.output_path + '/'
                              + self.output_header + '.bam')

    def clean_up_files(self):
        list = os.listdir(self.tmp_dir)
        for item in list:
            if self.output_header.split('.')[0] in item:
                if item.endswith(".fastq") \
                        or item.endswith(".sam")\
                        or item.endswith(".sai")\
                        or item.endswith("unsorted.bam"):
                    os.remove(join(self.tmp_dir, item))

    def run_alignment(
            self,
            no_index=False,
            tmp_dir='.',
            mark_dups=False,
    ):
        self.tmp_dir = tmp_dir
        self.logger.debug('Temp Dir:' + str(self.tmp_dir))
        self.logger.debug('No Index=' + str(no_index))
        self.logger.debug('Mark Duplicates=' + str(mark_dups))
        print "RUNNING ALIGNMENT..."
        # if build_index:
        #     self.build_index(bwa_mem=self.uses_mem())
        if not self.index_exists() and no_index:
            print "REFERENCE INDEX FILES REQUIRED BY BWA NOT PRESENT. ABORTING!"
            sys.exit(-1)
        if self.is_paired() and self.is_short():
            self.logger.debug('Run Alignment' + str(self.is_paired)
                              + ':' + str(self.is_short))
            if no_index is False:
                self.build_index(bwa_mem=False)
            self.__run_BWA_aln(1)
            self.__run_BWA_aln(2)
            self.__run_BWA_sam()
            self.__add_unmapped_reads(no_index)
            final_bam = self.tmp_dir + '/' + self.output_header + '.bam'
            if mark_dups:
                self.logger.info('Marking Duplicates')
                self.mark_duplicates(final_bam)
                final_bam = self.tmp_dir + '/' + self.output_header \
                            + '.duplicates_marked.bam'
            self.create_bam_index(final_bam)
        elif not self.is_paired() and self.is_short():
            self.logger.debug('Run Alignment' + str(self.is_paired)
                              + ':' + str(self.is_short))
            if no_index is False:
                self.build_index(bwa_mem=False)
            self.__run_BWA_aln(0)
            self.__run_BWA_sam()
            self.__add_unmapped_reads(no_index)
            final_bam = self.tmp_dir + '/' + self.output_header + '.bam'
            if mark_dups:
                self.logger.info('Marking Duplicates')
                self.mark_duplicates(final_bam)
                final_bam = self.tmp_dir + '/' + self.output_header \
                            + '.duplicates_marked.bam'
            self.create_bam_index(final_bam)
        elif not self.is_paired() and self.is_long():
            self.logger.debug('Run Alignment' + str(self.is_paired)
                              + ':' + str(self.is_long))
            if no_index is False:
                self.build_index(bwa_mem=False)
            self.__run_BWA_bwasw()
            self.__add_unmapped_reads(no_index)
            final_bam = self.tmp_dir + '/' + self.output_header + '.bam'
            self.create_bam_index(final_bam)
        elif self.uses_mem():
            self.logger.info('Using BWA MEM')
            if no_index is False:
                self.build_index(bwa_mem=True)
            self.__run_BWA_mem(paired=self.is_paired())
            self.__add_unmapped_reads(no_index, unsorted=True)
            final_bam = self.output_path + '/' + self.output_header \
                        + '.bam'
            aligned_bam = BamFile(self.tmp_dir + '/'
                                  + self.output_header + '_unsorted.bam'
                                  , 'bam', output=self.output_path + '/'
                                                  + self.output_header)
            aligned_bam.sort_bam()
            self.create_bam_index(final_bam)
        else:
            self.logger.info('Unable to run BWA on paired: '
                             + str(self.paired) + ' and short: '
                             + str(self.short) + ' Try using BowTie')
            sys.exit(-1)
        self.clean_up_files()


class BowTieAlignment(Alignment):

    def __init__(
            self,
            query,
            reference,
            output_header,
            paired=True,
            short=True,
            threads=2,
    ):
        full_path = os.path.dirname(os.path.abspath(output_header))
        full_header = os.path.basename(output_header)
        super(BowTieAlignment, self).__init__(
            query,
            reference,
            paired,
            short,
            threads,
            full_header,
            full_path,
        )
        self.index_name = os.path.basename(self.reference)

    def build_index(self):
        command(constant.BOWTIE + '/bowtie2-build ' + self.reference
                + ' ' + os.path.basename(self.reference))

    def build_bowtie_cmd(self, fastq_list, options=''):

        if len(fastq_list) == 1:
            read_string = '-U ' + fastq_list[0]
        elif len(fastq_list) == 2:
            read_string = '-1 ' + fastq_list[0] + ' -2 ' + fastq_list[1]
        else:
            print 'There are ' + str(len(fastq_list)) \
                  + 'fastq files in the list.  There should be 1 for an unpaired BAM or 2 for a paired BAM'
            sys.exit(-1)

        command = constant.BOWTIE + '/bowtie2 -p ' + str(self.threads) \
                  + ' ' + options + ' -x ' + self.index_name + ' ' \
                  + read_string + ' -S ' + self.output_path + '/' \
                  + self.output_header + '_bowtie.sam'

        return command

    def create_sorted_bam(self):
        aligned_sam = BamFile(self.tmp_dir + '/' + self.output_header
                              + '_bowtie.sam', 'sam', self.tmp_dir + '/'
                              + self.output_header + '_unsorted.bam')
        aligned_sam.convert_to_bam()

        aligned_bam = BamFile(self.tmp_dir + '/' + self.output_header
                              + '_unsorted.bam', 'bam',
                              self.output_path + '/'
                              + self.output_header)
        aligned_bam.sort_bam()

        self.create_bam_index(self.output_path + '/'
                              + self.output_header + '.bam')

    def clean_up_files(self):
        list = os.listdir(self.tmp_dir)
        for item in list:
            if self.output_header.split('.')[0] in item:
                if item.endswith(".fastq") \
                        or item.endswith(".sam") \
                        or item.endswith(".sai") \
                        or item.endswith("unsorted.bam"):
                    os.remove(join(self.tmp_dir, item))

    def run_alignment(self, no_index=False, tmp_dir='.'):

        self.tmp_dir = tmp_dir
        if no_index is False:
            self.build_index()
        fastq_list = self.convert_to_fastq(tmp_dir)
        if self.is_paired() and self.is_short():
            options = ''
        elif not self.is_paired() and self.is_short():
            options = ''
        elif not self.is_paired() and not self.is_short():
            options = '--local'
        elif self.is_paired() and not self.is_short():
            options = ''
        else:
            print 'Unable to run BowTie on paired: ' + str(self.paired) \
                  + ' and short: ' + str(self.short) + ' Try using BWA'
            sys.exit(-1)

        command(self.build_bowtie_cmd(fastq_list, options))
        self.create_sorted_bam()
        self.clean_up_files()
