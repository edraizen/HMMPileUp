#!/usr/bin/python
#Author: Eli Draizen
#Date: 4/25/2015
#File: run.py

#Standard libraries
import os, sys
import argparse
from datetime import datetime
from subprocess import Popen, PIPE

#Custom libraries
from HMMERPileUp import HMMERPileUp, is_hmmer_hmmsearch_file
from SAMPileUp import SAMPileUp, is_sam_output_directory

def run(files, reference_data, program="hmmer", filter_percentX=None, filter_spanX=None):
    if program == "hmmer":
        parser = HMMERPileUp(files[0])
    elif program == "sam":
        parser = SAMPileUp(files)
    else:
        raise RuntimeError("program must be 'hmmer' or 'sam'")
    parser.run(filter_percentX=filter_percentX, filter_spanX=filter_spanX)
    parser.printCysCount()
    parser.addGapsToRefData(reference_data)    
    parser.Write2File()
    return parser

def getHMMOutputFromDirectory(path, n=0):
    for f in os.listdir(path):
        if os.path.isdir(os.path.join(path,f)) and n<=2:
            openDirectory(os.path.join(path,f))
        elif f.endswith("hmmer_hmmer.out") or f.endswith("sam_hmmer.out"):
            yield [os.path.join(path,f)], "hmmer"
        elif f.endswith(".mstat") and is_sam_output_directory(os.path.join(path, f)[:-6]):
            yield is_sam_output_directory(os.path.join(path, f)[:-6], files=True), "sam"
        elif f.endswith(".mult"):
            continue
        elif is_hmmer_hmmsearch_file(os.path.join(path,f)):
            yield [os.path.join(path,f)], "hmmer"

def parse_args():
    """Parsing command line options
    """
    parser = argparse.ArgumentParser(description="")
    
    #Define Input
    parser.add_argument("reference_data",
                        help="Fasta file containing sequences or reference data")
    pipeline_options = parser.add_mutually_exclusive_group(required=True)
    pipeline_options.add_argument("--database",
                        default=None,
                        help="Path to database of fasta sequences for use in hmmsearch")
    pipeline_options.add_argument("--results",
                        default=None,
                        help="Path to hmm output if the hmm pipeline has already been run")
    parser.add_argument("--filter_percentX",
                        type=float,
                        required=False,
                        default=None,
                        help="Remove sequences that have have more then a given percentage of X's")
    parser.add_argument("--filter_spanX",
                        type=float,
                        required=False,
                        default=None,
                        help="Remove sequences that have a given number of consecutinve X's as percent")
 
    #Define output
    parser.add_argument("-o", "--outdir",
                        default="",
                        help="Path to save hmm output and pile ups. Default current directory")
    #Define output
    parser.add_argument("-l", "--log",
                        default="log.txt",
                        type=argparse.FileType('w'),
                        help="Path to save hmm output and pile ups. Default current directory")

    #Parse args
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    if args.database:
        raise RuntimeError("Error: cannot run pipline internally yet. Please run manually by editting the Makefile with correct parameters, and rerunning this script with the --results option.")

    print >> args.log, "Output from HMMPileUp.py"
    print >> args.log, "Run on", datetime.now()
    print >> args.log, "Arguments:", ";".join(["{}={}".format(name, value) for name, value in args.__dict__.iteritems()])
    
    if os.path.isdir(args.results):
        for files, program in getHMMOutputFromDirectory(args.results):
            result = run(files, args.reference_data, program=program, filter_percentX=args.filter_percentX, filter_spanX=args.filter_spanX)
            print >> args.log, program, files[0]
            print >> args.log, "  Total sequences:", result.total_seqs
            print >> args.log, "  Saved sequences:", len(result.records)
