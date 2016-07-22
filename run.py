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
from HMMPipeline import HMMPipeline
from HMMERPileUp import HMMERPileUp, is_hmmer_hmmsearch_file
from SAMPileUp import SAMPileUp, is_sam_output_directory

def run(files, reference_data, program="hmmer", filter_percentX=None, filter_spanX=None, filter_percentChar=None, filter_evalue=None):
    if program == "hmmer" and isinstance(files, [[], ()]):
        parser = HMMERPileUp(files[0])
    elif program == "hmmer" and isinstance(files, str):
        parser = HMMERPileUp(files)
    elif program == "sam" and isinstance(files, [[], ()]):
        parser = SAMPileUp(files)
    else:
        raise RuntimeError("program must be 'hmmer' (with 1 file) or 'sam' with (3 files)")
    parser.run(filter_percentX=filter_percentX, filter_spanX=filter_spanX, filter_percentChar=filter_percentChar, filter_evalue=filter_evalue)
    #parser.printCysCount()
    parser.addGapsToRefData(reference_data)    
    #parser.Write2File()
    return parser

def run_with_hmm(reference_data, database, filter_percentX=None, filter_spanX=None, filter_percentChar=None, filter_evalue=None):
    for (train, search), files in HMMPipeline.run(reference_data, database):
        result = run(
            files, 
            args.reference_data, 
            program=search, 
            filter_percentX=args.filter_percentX, 
            filter_spanX=args.filter_spanX, 
            filter_percentChar=args.filter_percentChar, 
            filter_evalue=args.filter_evalue
            )
        yield (train, search), result

def getHMMOutputFromDirectory(path, n=0):
    for f in os.listdir(path):
        if os.path.isdir(os.path.join(path,f)) and n<=2:
            for f in getHMMOutputFromDirectory(os.path.join(path,f)):
                yield f
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
    parser.add_argument("--filter_percentChar",
                        type=float,
                        required=False,
                        default=None,
                        help="Remove sequences that have have more then a given percentage of uppercase letters, not X")
    parser.add_argument("--filter_evalue",
                        type=float,
                        required=False,
                        default=None,
                        help="Remove sequences with evalue >= to value")
 
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
        

    print >> args.log, "Output from HMMPileUp.py"
    print >> args.log, "Run on", datetime.now()
    print >> args.log, "Arguments:"
    for name, value in args.__dict__.iteritems():
        print >> args.log, "  {}={}".format(name, value) 
    print >> args.log, ""
    
    if os.path.isdir(args.results):
        for files, program in getHMMOutputFromDirectory(args.results):
            result = run(
                files, 
                args.reference_data, 
                program=program, 
                filter_percentX=args.filter_percentX, 
                filter_spanX=args.filter_spanX, 
                filter_percentChar=args.filter_percentChar, 
                filter_evalue=args.filter_evalue
                )
            print >> args.log, program, files[0]
    elif args.database is not None:
        for (train, search), files in HMMPipeline.run(args.reference_data, args.database):
            result = run(
                files, 
                args.reference_data, 
                program=search, 
                filter_percentX=args.filter_percentX, 
                filter_spanX=args.filter_spanX, 
                filter_percentChar=args.filter_percentChar, 
                filter_evalue=args.filter_evalue
                )
            print >> args.log, search, files[0]
    print >> args.log, "  Total sequences:", result.total_seqs
    print >> args.log, "  Saved sequences:", len(result.records)
