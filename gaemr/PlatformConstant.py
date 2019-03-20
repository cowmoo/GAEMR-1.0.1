#!/usr/bin/env python

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


def get_bwa_version(path):
    proc = subprocess.Popen(path + "/bwa", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()
    err_split = err.split('\n')
    for item in err_split:
        if item.startswith("Version"):
            return item.split(":")[-1].rstrip("\n").lstrip(" ")
        # For some reason, when this code is uncommented this function returns None even when
        # when the condition is true. Makes no sense but debugging GAEMR is a hassle, so just commenting out for now.
        # else:
        #     return None


class PlatformConstant(object):
    """A unified class to contain all platform constants"""

    MIN_GAP_SIZE=10 # This constant is the minimum number of "N"s that constitute a gap when reading fasta files
    MIN_OUTPUT_GAP_SIZE=100  # This constant is in response to needing output gap sizes of 100 bases per NCBI
    QUAL_MAX=93
    MAX_INSERT_SIZE=50000
    TABLE_DELIMITER=" | "
    MUMMER_PATH="/broad/software/groups/gtba/software/mummer_3.23-64bit/" # "/path/to/mummer/package/" # http://mummer.sourceforge.net/
    BLAST_DIR="/broad/software/groups/gtba/software/ncbi-blast-2.2.25+/bin/" # "/path/to/blast+/package/" # http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download


    SAMTOOLS="/broad/software/groups/gtba/software/samtools_0.1.18/bin/samtools" # "/path/to/samtools/samtools" # http://samtools.sourceforge.net/
    #PICARD="/seq/software/picard/1.782/bin/" # "/path/to/public/picard/package/" # http://sourceforge.net/projects/picard/
    PICARD="/cil/shed/apps/external/picard/1.782/bin/"
    #PICARD="/seq/software/picard/current/bin/"
    #PICARD="/cil/shed/apps/external/picard/1.782/bin/" # "/path/to/public/picard/package/" # http://sourceforge.net/projects/picard/
    PICARD_TMP_DIR="."
    BLAST_NT="/broad/data/blastdb/nt/nt" # "nt"
    BLAST_UNIVEC="/cil/assembly_pipeline/databases/UniVec/UniVec" # "/path/to/UniVec/db/db" # ftp://ftp.ncbi.nih.gov/pub/UniVec/
    BLAST_rRNA="/cil/assembly_pipeline/databases/NCBI_rRNA/ncbi_rRNA" # "/path/to/rRNA/db/db" # Manually curated
    BLAST_MITOGCONTAM="/cil/assembly_pipeline/databases/mitogcontam/mitogcontam" # "/path/to/mitogcontam/db/db" # ftp://ftp.ncbi.nih.gov/pub/kitts/gcontam1.gz & ftp://ftp.ncbi.nih.gov/refseq/release/mitochondrion/*.f*a.gz
    BLAST_MITO="/cil/assembly_pipeline/databases/mitochondrion/mito_nt" # "/path/to/mito/db/db" # ftp://ftp.ncbi.nih.gov//blast/db/FASTA/mito.nt.gz
    # IMP: Changes this for production after we are happy with the new gcontam1
    # BLAST_GCONTAM="/cil/assembly_pipeline/databases/gcontam1/gcontam1" # "/path/to/mitogcontam/db/db" # ftp://ftp.ncbi.nih.gov/pub/kitts/gcontam1.gz
    BLAST_GCONTAM="/cil/assembly_pipeline/databases/gcontam1_staging/gcontam1" # "/path/to/mitogcontam/db/db" # ftp://ftp.ncbi.nih.gov/pub/kitts/gcontam1.gz
    BLAST_ADAPTORS="/cil/assembly_pipeline/databases/adaptors_for_screening/adaptors_for_screening" # NCBI's adaptors_for_screening database of next-generation sequencing primers
    BLAST_COMMON_PROK="/cil/assembly_pipeline/databases/common/prok/contam_in_prok.fa"
    BLAST_COMMON_EUK="/cil/assembly_pipeline/databases/common/euk/contam_in_euks.fa"
    BLAST_NODES="/broad/data/taxonomy/taxdump/nodes.dmp" # "/path/to/taxdump/nodes.dmp" # ftp://ftp.ncbi.nih.gov/pub/taxonomy/
    BLAST_NAMES="/broad/data/taxonomy/taxdump/names.dmp" # "/path/to/taxdump/names.dmp" # ftp://ftp.ncbi.nih.gov/pub/taxonomy/
    RNAMMER="/seq/annotation/bio_tools/rnammer/current/rnammer" # "/path/to/rnammer/package/rnammer" # http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer
    RDP="/broad/software/groups/gtba/software/rdp_classifier_2.4/rdp_classifier-2.4.jar" # "/path/to/rdp_classifier_2.4/package/rdp_classifier-2.4.jar" # http://sourceforge.net/projects/rdp-classifier/

    BLASTDBCMD=BLAST_DIR + "blastdbcmd"
    MAKEBLASTDB=BLAST_DIR + "makeblastdb"
    NUCMER=MUMMER_PATH + "nucmer"
    PROMER=MUMMER_PATH + "promer"
    SHOWTILING=MUMMER_PATH + "show-tiling"
    MUMMERPLOT=MUMMER_PATH + "mummerplot"

    BWA_MEM="/seq/software/picard/current/3rd_party/bwa_mem" # "/path/to/bwa/package/" # http://bio-bwa.sourceforge.net/
    BWA="/seq/software/picard/current/3rd_party/bwa"
    BWA_VERSION = get_bwa_version(BWA)
    BWA_MEM_VERSION = get_bwa_version(BWA_MEM)
    #MINIMAP2="/cil/shed/sandboxes/tshea/software/minimap2/minimap2"
    MINIMAP2="/cil/shed/apps/external/nanopore_software/general/bin/"
    GAEMR="/gsap/assembly_analysis/GAEMR/bin/" # /your/local/install/GAEMR
    BOWTIE="/broad/software/free/Linux/redhat_5_x86_64/pkgs/bowtie2_2.0.0-beta5/" # "/path/to/bowtie/package/" # http://bowtie-bio.sourceforge.net/index.shtml
    BOWTIE_VERSION='2.0.0-beta5'
    PLOT_COLORS=['green','blue','red','orange','magenta','black','grey']
    
    BLAST_XML_NO_DEF_LINE='No definition line'



    #detect user, environment and populate field
    def __init__(self):
        """populate class fields"""
        self.user_name = os.getenv("USER")
        self.operating_system = os.getenv("OSTYPE")


