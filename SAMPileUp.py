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

class SAMPileUp(HMMPileUp):
    def __init__(self, resultsFiles, reference_data, evalcutoff=0.0):
        HMMPileUp.__init__(self, resultsFiles[0], evalcutoff)
        self.mstatFile = resultsFiles[0] #Get HMM Length, evalues
        self.multFile = resultsFiles[1] #Get sequences
        self.resultsFiles = self.mstatFile
        self.hmmLength = len(SeqIO.parse(reference_data, "fasta").next())

    def parse(self, percentX=None):
        """
        """
        evalues = {}
        getEvalues = False
        with open(self.mstatFile) as mstat:
            for i, line in enumerate(mstat):
                if "% Sequence ID" in line:
                    getEvalues = True
                    continue
                if getEvalues:
                    sys.stdout.write("\r" + "."*(i%4))
                    sys.stdout.flush()
                    try:
                        _id, length, simple, _reversed, evalue = line.split()[:5]
                    except:
                        raise RuntimeError("Invalid mult file")

                    if float(evalue)<self.evalcutoff:
                        continue
                    evalues[_id] = float(evalue)

        if self.hmmLength == 0:
            raise RuntimeError("HMM Length cannot be be 0. Please check mstat file: {}".format(self.mstatFile))

        self.total_gaps = [0]*self.hmmLength

        #Assign evalues to sequences
        seqIDregex = re.compile("^(.+\|([0-9]+))_([0-9]+):([0-9]+)$")
        for i, seq in enumerate(SeqIO.parse(self.multFile, "fasta")):
            header_match = seqIDregex.search(seq.id)
            seqID = header_match.group(0)
            origSeqLength = header_match.group(1)
            seq = SAMSequence(seq.seq, self.hmmLength, origSeqLength)
            seq.align()
            seq.determineGapPositions()

            seqFrom = header_match.group(2)
            seqTo = header_match.group(3)
            seq.evalue = evalues[seqID]

            _id = "%d_%s" % (i+1, seqID)
            desc = "{}-{} [evalue: {}; program=UCSC-SAM3.5]".format(seqFrom, seqTo, seq.evalue)

            record = SeqRecord(seq, id=_id, description=desc)

            #Update gaps for all sequences, even if not saved
            self.updateGaps(seq.gaps)

            if not percentX or not seq.skip(percentX):
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

        self.data = MutableSeq("{}{}{}".format(begin, matchedSeq, end)).data
