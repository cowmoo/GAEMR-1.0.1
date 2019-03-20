#
# Class for manipulating bam files
#

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
import os
from RunCommand import RunCommand
from FastqFile import FastqFile
import PlatformConstant as pc
import subprocess


constant = pc.PlatformConstant()

class BamFile():
    """ This class represents a manipulation of BAM files. """

    def __init__(self, file, extension, output=None, direction=None, sample='sample'):
        self.file = file
        self.extension = extension
        self.direction = None
        if output:
            self.output = output
        else:
            self.output = self.__get_output_filename()

        if direction:
            self.direction = direction.upper()
        
        self.sample = sample
            

    #converts BAM file reads into a Fasta file of reads
    def convert_to_fasta(self,fasta_file_prefix):
        """converts BAM file reads into a Fasta file of reads"""
        fasta_conversion_cmd = ['java','-jar', constant.PICARD + "picard.jar", "SamFormatConverter"]
        fasta_conversion_cmd.append("INPUT=" + self.file)
        samfile_name  = self.file.replace('.bam','.sam')
        fasta_conversion_cmd.append("OUTPUT=" + samfile_name)
        rc = RunCommand(fasta_conversion_cmd)
	print fasta_conversion_cmd
	out = rc.run_command()
	print out
	
        #convert sam file into fasta file of reads
        sam_file = open(samfile_name,'r')
        fasta_file_name = fasta_file_prefix + '.fasta'
        fasta_file = open(fasta_file_name,'w')
        for line in sam_file:
	    if line[0] != '@':
            	fields = line.split('\t')
            	name = fields[0]
            	seq = fields[9]
            	fasta_file.write('>')
            	fasta_file.write(name)
            	fasta_file.write('\n')
            	fasta_file.write(seq)
            	fasta_file.write('\n')
        fasta_file.close()
    	return fasta_file_name
                
                    

    def __get_output_filename(self):
        return str(self.__get_base_header()) + ".unmapped.bam"
    
    def __get_base_header(self):
        names = self.get_input_file().split('.')
        return names[0] + "." + names[1]

    def __get_direction(self):
        return self.direction

    def get_output_file(self):
        return self.output

    def __build_fastq_command(self, java_heap='32G', stringency='LENIENT'):
        output = self.get_output_file()
        tmp_dir = "TMP_DIR=" + constant.PICARD_TMP_DIR
        fastq_command = constant.PICARD + "picard.jar"

        read1_fq = re.sub("bam","read1.fastq",output)
        read2_fq = re.sub("bam","read2.fastq",output)
        fastqs = [read1_fq,read2_fq]
        fastq_list = ["F=",read1_fq,"F2=",read2_fq]
        if not self.__get_direction():
            fastqs = [read1_fq]
            fastq_list = ["F=",read1_fq]
        return ["java","-Xmx" + java_heap,"-jar",fastq_command,"SamToFastq", "NON_PF=true",tmp_dir,
                "I=",self.get_input_file(), "VALIDATION_STRINGENCY=", stringency] + fastq_list,fastqs

    def __remove_fastq_files(self,fastqs):
        for i in fastqs:
            os.remove(i)

    def convert_to_fastq(self, java_heap='32G', stringency='LENIENT'):
        cmd,fastqs = self.__build_fastq_command(java_heap, stringency)
        if cmd:
            rc = RunCommand(cmd)
            out = rc.run_command()
        return 1,fastqs

    def convert_to_unmapped_bam(self, java_heap='32G', stringency='LENIENT'):

        direction = self.__get_direction()
        if not direction or direction == 'FR':
            cmd = self.__build_command(java_heap, stringency)
            if cmd:
                rc = RunCommand(cmd)
                out = rc.run_command()
            else:
                if self.is_mapped and not os.path.exists(self.get_output_file()) and self.get_output_file() != self.get_input_file():
                    os.symlink(self.get_input_file(),self.get_output_file())
            return 1
        else:
            fq_return_code, fastqs = self.convert_to_fastq(java_heap, stringency)
            fq = FastqFile(fastqs,self.sample,None,self.__get_direction(),self.get_output_file())
            fq_return_code = fq.convert_to_unmapped_bam(java_heap, stringency)
            self.__remove_fastq_files(fastqs)
            return fq_return_code

    def __build_command(self,java_heap='32G', stringency="LENIENT"):
        convert_command = constant.PICARD + "picard.jar"
        tmp_dir = "TMP_DIR=" + constant.PICARD_TMP_DIR
        if self.is_mapped():
            return ["java","-Xmx" + java_heap,"-jar",convert_command, "RevertSam", "I=", self.get_input_file(),
                    "O=", self.get_output_file(), "VALIDATION_STRINGENCY=", stringency, 'SAMPLE_ALIAS=', self.sample,
                    "OQ=true", "REMOVE_DUPLICATE_INFORMATION=true", "REMOVE_ALIGNMENT_INFORMATION=true",tmp_dir]
        return 0
    
    def __check_headers(self,out):
        mapped = 0
        for i in out.split('\n'):
            if re.match("^@SQ",i):
                mapped = 1
        return mapped
    
    def is_mapped(self):
        cmd = None
        if self.__is_sam_format():
            cmd = ["grep","^@",self.get_input_file()]
        elif self.__is_bam_format():
            cmd = [constant.SAMTOOLS,"view", "-H",self.get_input_file()]
        else:
            return 0
        rc = RunCommand(cmd)
        out = rc.run_command()
        return self.__check_headers(out)

    def __is_sam_format(self):
        type = self.__get_input_filetype()
        return re.match("SAM",type) or re.match("Sam",type) or re.match("sam",type)
    
    def __is_bam_format(self):
        type = self.__get_input_filetype()
        return re.match("BAM",type) or re.match("Bam",type) or re.match("bam",type)

    def get_input_file(self):
        return self.file

    def __get_input_filetype(self):
        return self.extension

    def convert_to_bam(self):
        if self.__is_bam_format():
            print self.get_input_file() + " is already a BAM file"
            return 0
        else:
            #s = "-S"
            #cmd = [constant.SAMTOOLS, "view", "-b", s, self.get_input_file(),">", self.get_output_file()]
            cmd = constant.SAMTOOLS + " view -b -S " + self.get_input_file() + " > " + self.get_output_file()
            self.command(cmd)
            #rc = RunCommand(cmd)
            #print rc.get_command()
            #out = rc.run_command()
        return 1

    def sort_bam(self):

        cmd = [constant.SAMTOOLS, "sort", self.get_input_file(), self.get_output_file()]
        rc = RunCommand(cmd)
        out = rc.run_command()
        return 1


    def command(self,cmd):
        print cmd
        subprocess.check_call(cmd, shell=True)
