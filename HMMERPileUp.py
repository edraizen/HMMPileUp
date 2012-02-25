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

# Modules from Biopython
from Bio.Seq import Seq, MutableSeq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import Bio.HMMER as HMM #Must have SPA's biopython fork

# Local modules
from Sequences import Sequences, MutableSequences

# Part One -- Parsing the hmmsearch file
 
class HMMERPileUp(object):
    """These are a collection of scripts that convert the output from hmmer3's
    hmmsearch to aligned FASTA with the included gaps. It also adds relevant
    info of the hit into the FASTA header and optionally adds all of the gaps
    into the reference data."""
    
    def __init__(self, resultsFile):
        """Initialize the parser.

        Parameters:
        resultsFile - string. path/to/hmmsearch_output.txt
        """
        self.resultsFile = resultsFile
        self.hmmLength = 0
        self.records = []
        self.total_gaps = []
        self.bitThreshold = 0.0
        self.refDB = None
        
    def parse(self, bitThreshold = 0, align=True, gaps=True):
        """Search through all of the hits in the given hmmer3 output.
        Hits are stored into an HMMSequence object, which counts
        the gaps in each hit. All of the positions and lengths of the gaps are  stored into the
        parser to add into all of the hits. Information about the hit is
        also saved, and will be turned into a FASTA header.

        Parameters:
        bitThreshold - a floating point value. Ignore bits that are less than bitThreshold.
        align - boolean. Enables alignment algorithm based on HMM values.
                See HMMSequence.align() for more info.
        gaps - boolean. Enables gap counting algorithm to create pileup.
               See HMMSequence.determineGapPositions()
        """
        
        hmmFile = open(self.resultsFile)
        parsedHMM = HMM.parseHMMER3(hmmFile)
        self.hmmLength = int(parsedHMM.hmmLength)
        self.total_gaps = [0]*self.hmmLength
        self.bitThreshold = bitThreshold
        
        removedSequences = []
        
        for unit in parsedHMM:
            bits = unit.bits
            if bits < bitThreshold:
                continue

            seq = HMMSequence(unit, self.hmmLength, align=align, gaps=gaps)
            _id = "%d_%s" % (len(self.records)+1, unit.name)
            desc = "%d-%d (HMM: %d-%d) [score: %0.1f bits]" % (unit.seqFrom,
                                                              unit.seqTo,
                                                              unit.hmmFrom,
                                                              unit.hmmTo,
                                                              unit.bits)
            record = SeqRecord(seq, id=_id, description=desc)
            
            #Remove sequence if there is more 50% X in the string
            if record.seq.count("X")>len(record.seq)/2:
				removedSequences.append(record.description)
				continue
				
            self.records.append(record)
            self.__addGaps(seq.gaps)
        
        if len(removedSequences) > 0:
            print "Removed sequences from %s:" % os.path.basename(self.resultsFile)
            for seq in removedSequences:
                print seq
        
    def __addGaps(self, gaps):
        """Private method that adds the gaps from an HMMSequence
        to this parser. Gaps are stored in an index with index as
        the postion of the sequence to add the gaps and the value
        as the number of gaps at that postion. The array is currently
        the length of the hmm and if there is no gap at a given index,
        the value at that index is zero.
        """
        
        for index, gap in enumerate(gaps):
            if gap > self.total_gaps[index]:
                self.total_gaps[index] = gap
                
    def setRefDB(self, file):
    	self.refDB = file

    def addGapsToRefData(self, name=None, format="fasta"):
        """This method is where all of the gaps are inserted into
        the reference data.

        Parameters:
        refData - string. path/to/refData.fasta. Accepts FASTA and other
                  Biopython compatible formats.
        name - string. Name to save the ref data, if not same name. Optional
        format - string. Format of refData. Accepts FASTA and other
               Biopython compatible formats. Optional if format is FASTA.
        """
        
        if name is None:
        	name = "RefData_%s.txt" % (os.path.basename(self.resultsFile))
        #if self.bitThreshold>0.0:
        #name = "%s_bitsThreshold_%d" % (".", self.bitThreshold)
        print self.refDB    
        sequences = MutableSequences(self.refDB, format=format)
        realGaps = self.total_gaps[:]
        realGaps.reverse()
        for i, seq in enumerate(sequences):
            seq.id = "%d_%s" % (i+1, seq.id)
            for index, gap in enumerate(realGaps):
                index = self.hmmLength-index-1
                if gap == 0: continue
                for i in range(gap):
                    seq.seq.insert(index+1, "-")
        SeqIO.write(sequences, name, format)
        
    def addGapsToHMMSeqs(self):
        """Re-inserts all gaps into hit sequences. For more info
        see HMMSequence.insertAllGaps()
        """
        
        for seq in self.records:
            seq.seq.insertAllGaps(self.total_gaps)
        
    def ignoreBitsBelow(self, bitThreshhold):
        """Ignore bits below a certain threshhold. Not needed if supplied threshold as an input.

        Parameters:
        bitThreshold - a floating point value. Ignore bits that are less than bitThreshold.
        """
        for record in records:
            if record.hmm.bits < threshold:
                del record

    def Write2File(self, name=None, type="fasta"):
        """Save the hit sequences to any format

        name - string. Name to save the ref data, if not same name.
        type - string. Format of refData. Accepts FASTA and other
               Biopython compatible formats. Optional if format is FASTA.
        """
        
        if name is None:
            name = "%s_PILEUP.txt" % (os.path.splitext(os.path.basename(self.resultsFile))[0])
            
        tmp = []
        for i, seq in enumerate(self.records):
            tmp_seq = seq
            if not seq.__class__.__name__ == "SeqRecord":
                tmp_seq.seq = seq.seq.toseq()
            tmp.append(tmp_seq)
        SeqIO.write(tmp, name, type)

    def WriteGaps2File(self, name=None):
        """Save the global gap list in a tab delimeted file

        Parameters:
        name - string. Name to save the gap list. Optional.
               If not included, name with the name of the results file ending wiht with "GAPS.txt"
        """

        if name is None:
            name = "%s_GAP_List.txt" % os.path.basename(self.resultsFile).split(".")[0]
            
        f = open(name, 'w')
        f.write("index\tlength")
        for index, gap in enumerate(self.total_gaps):
            if gap == 0: continue
            f.write(str(index) + "\t" + str(gap))
        f.close()

    def printTotalGaps(self):
        """Prints the global gap list to the standard out
        """
        
        print "index\tgap length"
        for index, gap in enumerate(self.total_gaps):
            if gap == 0: continue
            print str(index) + "\t" + str(gap)

# Part Two -- Align the Hit sequence to the HMM model 
# and Insert Gap positions from the HMM model

# For unknown amino acids at each position not accounted for in the HMM model
# (SPA is 124 positions per domain) insert an "X". 

# For each gap position in the HMM model (indicated by a .) insert a
# corresponding "-". As gap lengths at each position can vary, insert the
# HIGHEST gap length value for that position.
 
class HMMSequence(MutableSeq):
    """An object that bridges an HMM hit sequence with Biopython's Seq object.
    This object determines the gaps of the hmm and insert them into the sequence
    and saves all of the information from the HMM hit into the sequence."""
    
    def __init__(self, unit, hmmLength, align=True, gaps=True):
        """Intialise HMMSequence with the hmmer unit.

        Parameters:
        unit - HMMUnit object. 
        hmmLength - int. length of the HMM.
        align - boolean. Enables alignment algorithm based on HMM values.
                See HMMSequence.align() for more info.
        gaps - boolean. Enables gap counting algorithm to create pileup.
               See HMMSequence.determineGapPositions()
        """
        
        super(HMMSequence, self).__init__(unit.hmmalign["seq"])
        self.hmm = unit
        self.hmmLength = int(hmmLength)
        self.gaps = [0]*self.hmmLength
        if align:
            self.align()
        if gaps:
            self.determineGapPositions()

    def addToFront(self, seq):
        """Insert a sequence in front of this sequence.

        Parameters:
        seq - A Sequence Object or string. The sequence to add.
        """
        
        for c in reversed(str(seq)):
            self.data.insert(0, c)
            
    def insertUnknowns(self, index, unknowns):
        """Insert a certain number of unknows into this sequence, denoted by an "X"

        Parameters:
        index - int. The position to insert the unknowns in teh sequence.
        unknonws  - int. The number of unknonws to add into the sequence. 
        """
        
        for i in range(unknowns):
            self.data.insert(index, "X")

    def insertDeletionsAfter(self, index, deletions):
        # Here if we encounter a lower case letter we
        # subtract 1 from the deletions

        #i = prevEnd
        #while i >= index:
        #    AtIndex = self.data[i]
        #    if (not (AtIndex == '-' or AtIndex.isupper())):
        #        deletions-=1
        #    i-=1
        #if deletions > 0:
        for i, j in enumerate(range(deletions)):
            self.data.insert(index+1, "-")

    def align(self):
        """Add gaps and unknowns to this HMM sequence based
        on information of this HMM hit
        """
        number_of_Xs = 0
        xFront = ""
        xEnd = ""
        dashFront = ""
        dashEnd = ""
        domain = self.hmm.domain
        seqFrom = self.hmm.seqFrom
        seqTo = self.hmm.seqTo
        hmmFrom = self.hmm.hmmFrom
        hmmTo = self.hmm.hmmTo

        # Determining if variable amino acids ("X") need to be added to the
	    # beginning of the sequence:
        z = hmmFrom-seqFrom
        number_of_Xs = (hmmFrom-1)-z
        if z > 0:
            dashFront = "-"*z
            xFront = "X"*number_of_Xs
        elif hmmFrom-1<seqFrom-1:
            xFront = "X"*(hmmFrom-1) 

        # Determining if variable amino acids ("X") need to be added to the 
        # end of the sequence:
        seqRange = seqTo-seqFrom
        number_of_Xs_end = self.hmmLength - hmmTo

        # The original sequence length; SPA format includes this
        lengthPattern = re.compile(r'(\d+)$')
        res = lengthPattern.search(self.hmm.name)
        seqLength = int(res.group(1)) #int(self.hmm.name.split("_")[-1])
        distToSeqEnd = seqLength - seqTo
        if distToSeqEnd >= number_of_Xs_end and number_of_Xs_end != self.hmmLength:
            xEnd = 'X'*number_of_Xs_end
        else:
            if distToSeqEnd < number_of_Xs_end:
                xEnd = 'X'*distToSeqEnd
            dashEnd += "-"*(number_of_Xs_end-distToSeqEnd)

        begin = "%s%s" % (dashFront, xFront)
        end = "%s%s" % (xEnd, dashEnd)
        self.addToFront(begin)
        self.data.extend(end)
        self.original = str(self)

    def determineGapPositions(self):
        """Determine how many gap positions and the indexes
        of the gaps in this hit
        """
        
        total_letters = -1
        dots = 0
        total_dots = 0
        for c in self.hmm.hmmalign["hmm"]:
            if c != ".":
                total_letters += 1
            else:
                dots += 1
            if not c == "." and dots > 0:
                gap_index = total_letters+self.hmm.hmmFrom-2
                gap = dots
                self.gaps[gap_index] = gap
                total_dots += dots
                dots = 0

    def printGaps(self, gaps):
        print "index\tgap length"
        for index, gap in enumerate(gaps):
            if gap == 0: continue
            print str(index) + "\t" + str(gap)

    def insertAllGaps(self, gaps):
        """This method inserts the previously computed gaps starting at the
        end of the sequence by reversing the list and recomputing.
        """

        total_lower = sum(self.gaps)
        refIndex = self.hmmLength-1
        index = len(self.data)-1 #self.hmmLength-1
        while index >=0:
            refIndex = index-(total_lower-sum(self.gaps[refIndex:]))
            index -= self.gaps[refIndex]
            gap = gaps[refIndex]-self.gaps[refIndex]
            if gap>0:
                self.insertDeletionsAfter(index, gap)
            index-=1

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["bitThreshold=", "refData=", "onlyRef"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    except e:
    	print e
        
    bitThreshold = 0
    refData = None
    onlyRef = False
    
    for opt, arg in opts:
        if opt in ["--bitThreshold"]:
            bitThreshold = float(arg)
        elif opt in ["--refData"]:
            refData = arg
        elif opt in ["--onlyRef"]:
        	onlyRef = True
        else:
            assert False, "unhandled option"
    try:
        hmmFile = args[0]
    except:
        print "Invalid input"
        sys.exit(2)
    
    parser = HMMERPileUp(hmmFile)
    parser.parse(bitThreshold=bitThreshold)
    if refData is not None:
    	parser.setRefDB(refData)
    	parser.addGapsToRefData()
    if not onlyRef:
    	parser.addGapsToHMMSeqs()
    parser.Write2File()
    
    
