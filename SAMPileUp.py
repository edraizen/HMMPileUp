#!/usr/bin/python

#SAMPileUp.py by Eli Draizen (edraizen@soe.ucsc.edu)

## Expected input is an hmmscore file (in .mstat format) and a prettyalign  ##
## (in .mult format) command line argument.  Output is a file containing   ##
## FASTA formatted sequences from the hits that matched our HMMs including ##
## insertions for variable amino acids (as X's) and gap positions (as -'s) ##

# Modules from python standard
import os, sys
import re
import getopt

# Modules from Biopython
from Bio.Seq import Seq, MutableSeq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Local modules
from HMMPileUp import HMMPileUp, HMMSequence

def is_sam_output_directory(path, files=False):
    needed_files = ["{}.{}".format(path, ext) for ext in ("dist", "mstat", "mult")]

    if files:
        return needed_files

    return all(os.path.isfile(f) for f in needed_files)

class SAMPileUp(HMMPileUp):
    def __init__(self, resultsFiles, evalcutoff=0.0):
        HMMPileUp.__init__(self, resultsFiles[0], evalcutoff)
        self.distFile = resultsFiles[0] #Get HMM Length
        self.mstatFile = resultsFiles[1] #Get e-values
        self.multFile = resultsFiles[2] #Get sequences
        self.resultsFiles = self.mstatFile

    def get_hmm_length(self):
        """Get HMM Length fro Dist file.
        Match lines like:
        % 264051 sequences, 157742148 residues, 124 nodes, 1834.71 seconds
        """
        nodes_re = re.compile(r', (\d+) nodes,')
        with open(self.distFile) as dist:
            for line in dist:
                match = nodes_re.search(line)
                if match:
                    self.hmmLength = int(match.groups()[0])-1
                    return
        if self.hmmLength == 0:
            raise RuntimeError("Invalid SAM HMM, with length 0. Please check {} is accurate.".format(self.distFile))

    def parse(self):
        """
        """
        self.get_hmm_length()

        evalues = {}
        getEvalues = False
        program = "SAM"
        with open(self.mstatFile) as mstat:
            for i, line in enumerate(mstat):
                if line.startswith("% SAM"):
                    program = line[1:].rstrip()
                if line.startswith("% Sequence ID"):
                    getEvalues = True
                    continue
                if getEvalues:
                    try:
                        _id, length, simple, _reversed, evalue = line.split()[:5]
                    except:
                        raise RuntimeError("Invalid mult file")
                    if float(evalue)<self.evalcutoff:
                        continue
                    evalues[_id] = float(evalue)

        self.total_gaps = [0]*self.hmmLength

        #Assign evalues to sequences
        seqIDregex = re.compile("^(.+\|([0-9]+))_([0-9]+):([0-9]+)$")
        for i, s in enumerate(SeqIO.parse(self.multFile, "fasta")):
            header_match = seqIDregex.search(s.id)
            seqID = header_match.groups()[0]
            origSeqLength = header_match.groups()[1]

            seqFrom = header_match.groups()[2]
            seqTo = header_match.groups()[3]

            seq = SAMSequence(s.seq, self.hmmLength, origSeqLength, evalues[s.id], seqFrom, seqTo)
            seq.align()
            seq.determineGapPositions()

            

            _id = "%d_%s" % (i+1, seqID)
            _id = "{}_{}".format(i+1, seqID)
            desc = "[Seq:{}-{}; HMM: {}-{}; e-value: {}; program={}]".format(
                    seqFrom,
                    seqTo,
                    seq.hmm_start,
                    seq.hmm_end,
                    seq.evalue,
                    program
                    )

            record = SeqRecord(seq, id=_id, description=desc)

            #Update gaps for all sequences, even if not saved
            self.updateGaps(seq.gaps)

            if not seq.skip():
                self.records.append(record)

class SAMSequence(HMMSequence):
    seq_regex = re.compile('([a-z]*)(-*)([A-Z][A-Za-z\-]*[A-Z])(-*)([a-z]*)')
    def align(self):
        #extract dashes and capital letters
        #get rid of everything before upper and after
        #lowercase in front number
        #if lower case in middle it's an insertion
        #if dash
        #coutn dashes at start before uppercase
        #count lower case before UPPER or DASH
        #           0      1       2(count 1st)      3    4
        #reg ex: ([a-z]*)(-*)([A-Z][A-Za-z\-]*[A-Z])(-*)([a-z]*)
        # 1,2,3 count number
        #[A-Z] in front and back make sure there are upper case 
        m = SAMSequence.seq_regex.search(str(self))
        frontLowers = m.groups(0)[0]
        frontDashes = m.groups(1)[1]
        matchedSeq = m.groups(1)[2]
        endDashes = m.groups(1)[3]
        endLowers = m.groups(1)[4]
        
        #Append X's to beginning
        if len(frontDashes)<=len(frontLowers):
            begin = "X"*len(frontDashes)
        else:
            begin = "-"*(len(frontDashes)-len(frontLowers))
            begin += "X"*len(frontLowers)

        #Append X's to end
        if len(endDashes)<=len(endLowers):
            end = "X"*len(endDashes)
        else:
            end = "X"*len(endLowers)
            end += "-"*(len(endDashes)-len(endLowers))

        self.hmmStart = len(frontLowers)
        self.hmmEnd = len(self)-len(endLowers)+1

        self.data = MutableSeq("{}{}{}".format(begin, matchedSeq, end)).data
