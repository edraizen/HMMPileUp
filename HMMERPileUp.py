#!/usr/bin/python

# HMMERPileUp.py originally by Eli Draizen (edraizen@soe.ucsc.edu)
# contributions by Felicia Kemp (felicia@soe.ucsc.edu) and Eric Christopher (echristo@gmail.com)

## Expected input is an hmmsearch file (in .txt or .hmm format) given as a ##
## command line argument.  Output is a file containing FASTA formatted     ##
## sequences from the hits that matched our HMMs including insertions for  ##
## variable amino acids (as X's) and gap positions (as -'s).               ##


# Modules from python standard
import sys, os
import getopt
import re
import math
import warnings

# Modules from Biopython
from Bio import BiopythonExperimentalWarning
warnings.simplefilter('ignore', BiopythonExperimentalWarning)
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO

# Local modules
from HMMPileUp import HMMPileUp, HMMSequence

def is_hmmer_hmmsearch_file(f):
    """Return true if file is of an hmmsearch run"""
    #HMMER output doesn't have a specific file extension
    with open(f) as hmmer_file:
        try:
            return hmmer_file.next().startswith("% hmmsearch")
        except StopIteration:
            return False

# Part One -- Parsing the hmmsearch file
class HMMERPileUp(HMMPileUp):
    """These are a collection of scripts that convert the output from hmmer3's
    hmmsearch to aligned FASTA with the included gaps. It also adds relevant
    info of the hit into the FASTA header and optionally adds all of the gaps
    into the reference data."""
        
    def parse(self):
        """Search through all of the hits in the given hmmer3 output.
        Hits are stored into an HMMSequence object, which counts
        the gaps in each hit. All of the positions and lengths of the gaps 
        are  stored into the parser to add into all of the hits. Information
        about the hit is also saved, and will be turned into a FASTA header.
        """

        try:
            query = SearchIO.parse(self.resultsFile, "hmmer3-text").next()
        except StopIteration:
            raise RuntimeError("Invalid HMMER output")


        self.hmmLength = query.seq_len
        self.total_gaps = [0]*self.hmmLength
        num_hits = 0
        for i, hit in enumerate(query):
            #if not hit.is_included:
                #Skip sequences below threshold
                #continue
            origSeqLength = int(hit.id.split("|")[-1])
            for j, hsp in enumerate(hit):
                num_hits += 1
                seq = HMMERSequence(
                    str(hsp.hit.seq), 
                    query.seq_len, 
                    origSeqLength, 
                    hsp.evalue,
                    hsp.hit_start, 
                    hsp.hit_end, 
                    hsp.query_start, 
                    hsp.query_end
                    )
                seq.align(hsp.hit_start, hsp.hit_end, hsp.query_start, hsp.query_end)
                seq.determineGapPositions()
                _id = "{}_{}".format(num_hits, hit.id)
                desc = "[Seq:{}-{}; HMM: {}-{}; e-value: {}; program={}]".format(
                    hsp.hit_start+1,
                    hsp.hit_end,
                    hsp.query_start,
                    hsp.query_end,
                    hsp.evalue,
                    query.program
                    )
                record = SeqRecord(seq, id=_id, description=desc)

                #Update gaps for all sequences, even if not saved
                self.updateGaps(seq.gaps)

                if not seq.skip() and hit.is_included:
                    self.records.append(record)

# Part Two -- Align the Hit sequence to the HMM model 
# and Insert Gap positions from the HMM model

# For unknown amino acids at each position not accounted for in the HMM model
# (SPA is 124 positions per domain) insert an "X". 

# For each gap position in the HMM model (indicated by a .) insert a
# corresponding "-". As gap lengths at each position can vary, insert the
# HIGHEST gap length value for that position.
 
class HMMERSequence(HMMSequence):
    def align(self):
        """Add gaps and unknowns to this HMM sequence based
        on information of this HMM hit
        """
        number_of_Xs = 0
        xFront = ""
        xEnd = ""
        dashFront = ""
        dashEnd = ""

        # Determining if variable amino acids ("X") need to be added to the
	    # beginning of the sequence:
        z = self.hmmStart-self.seqStart
        number_of_Xs = (self.hmmStart-1)-z
        if z > 0:
            dashFront = "-"*z
            xFront = "X"*number_of_Xs
        elif self.hmmStart-1<=self.seqStart-1:
            xFront = "X"*(self.hmmStart-1) 

        # Determining if variable amino acids ("X") need to be added to the 
        # end of the sequence:
        number_of_Xs_end = self.hmmLength - self.hmmEnd

        # The original sequence length; SPA format includes this
        delimeter = "|" #Need to fix can be "_" or "|" or something else...
        
        distToSeqEnd = self.origSeqLength - seqTo
        if distToSeqEnd >= number_of_Xs_end and number_of_Xs_end != self.hmmLength:
            xEnd = 'X'*number_of_Xs_end
        else:
            if distToSeqEnd < number_of_Xs_end:
                xEnd = 'X'*distToSeqEnd
        	dashEnd += "-"*(number_of_Xs_end-distToSeqEnd)
        	
        begin = "{}{}".format(dashFront, xFront)
        end = "{}{}".format(xEnd, dashEnd)
        self.addToFront(begin)
        self.data.extend(end)
        self.original = str(self)
