"""
===========================
pipeline_rMATS
===========================

:Author: Jack Riley

Overview
========

Pipeline to detect differential alternative splicing from RNA-sequencing data
by using the rMATS software (https://github.com/Xinglab/rmats-turbo/blob/v4.1.1). 

Currently set up to run against a fixed event set (3'UTR introns). Although this
will later be made configurable to allow complete comparison. 

Usage
=====

To execute this pipeline in an interactive ShARC session or from a personal
device (not recommended) use the following command:

"python <path_to_pipeline_folder>/pipeline_rMATS.py make full -v5"

To execute this pipeline to the ShARC cluster, use the following:

"submit_pipeline <path_to_pipeline_folder>/pipeline_rMATS.py make full -v5"


Configuration
=============

In order for pipeline_rMATS to know what comparison to make, 2 important files are required
(without which the pipeline will not work): 

    - design.tsv = tsv file where each line a comparison of 2 named variables 
        - example design.tsv can be found at pipeline_rMATS/pipeline_rMATS/design.tsv
    - file_naming.tsv = tsv file which links the named variables to regex to the files of interest by
        - example file_naming.tsv can be found at pipeline_rMATS/pipeline_rMATS/file_naming.tsv


Input files
===========

Files are given as .bam format. Usually following completion of pipeline_utrons. (.fastq files can be
used with rMATS however this has not yet been configured). In order to compare against fixed event sets, 
files from utron_beds (output from pipeline_utrons) will also be needed.


Requirements
============

Packages required:
    - cgat/cgatcore
    - rmats-turbo (v4.1.1 was used during creation of pipeline)
        - a compatible conda environment has been created and
        stored in sudlab shared drive ('rmats-env')

These packages are installed in the conda environment "qapa-env".
R packages for final analysis/reports are installed in "qapa-env-R". 


Pipeline output
===============

rMATS produces tabular output for each type of alternative splicing event. These will be used by a custom
Rmarkdown to produce a human-readable summary of the rMATS results. Subsequent data manulation can then be
conducted manually.

Code
====

"""

###################
##### imports #####
###################

from ruffus import *
from ruffus.combinatorics import product
import sys
import os
import shutil
from gffutils import DataIterator as DataIterator
import sqlite3
import subprocess
import glob
import csv
from cgatcore import experiment as E
import cgat.Sra as Sra
from cgatcore import pipeline as P
import cgatpipelines.tasks.rnaseq as RnaSeq
import tempfile


############################
#####       main       #####
############################

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    'genesets',
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

PARAMS["project_src"]=os.path.dirname(__file__)

RnaSeq.PARAMS = PARAMS


###################
##### utility #####
###################

@follows(downloadQAPAprereqs, build3UTRlib, extract3UTRseq, makeSalmonIndex, 
        quantifyWithSalmon, quant3UTRusage, analyseDesignMatrix,
        run_compareQAPA_script)
def full():
    pass


##################
###### misc ######
##################

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))