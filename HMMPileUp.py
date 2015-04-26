#Standard libraires
import os
from math import ceil, floor

#Required libraires
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

#Custom libraires
from Sequences import MutableSequences

class HMMPileUp(object):
    def __init__(self, resultsFile, evalcutoff=0.0):
        self.resultsFile = resultsFile
        self.hmmLength = 0
        self.evalcutoff = evalcutoff
        self.total_gaps = [0]*self.hmmLength
        self.records = []

    def run(self, percentX=None):
        """Run the HMM parser and insert X's and dashes based on the 
        alignment and then adds the total gaps found in alignments to 
        all of the sequences. Then it removes sequences whose sequences
        contain less than percentAA number of amino acids.

        Parameters:
        percentAA - how many amino acids are allowed in the sequence
            default is 0.5
        """
        self.parse(percentX=percentX)
        self.addGapsToHMMSeqs()

    def parse(self):
        """
        """
        raise NotImplementedError

    def updateGaps(self, gaps):
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

    def addGapsToHMMSeqs(self):
        """Re-inserts all gaps into hit sequences. For more info
        see HMMSequence.insertAllGaps()
        """
        for seq in self.records:
            seq.seq.insertAllGaps(self.total_gaps)

    def addGapsToRefData(self, refDB, name=None, file_format="fasta"):
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
        	name = "RefData_{}_PILEUP.txt".format(os.path.splitext(os.path.basename(self.resultsFile))[0])
        
        reference_seqs = []
        for s in SeqIO.parse(refDB, file_format):
            s.seq = HMMSequence(str(s.seq), self.hmmLength, self.hmmLength)
            s.seq.insertAllGaps(self.total_gaps)
            reference_seqs.append(s)

        self.Write2File(sequences=reference_seqs, name=name, file_format=file_format)

    def Write2File(self, name=None, file_format="fasta", sequences=None):
        """Save the hit sequences to any format

        name - string. Name to save the ref data, if not same name.
        format - string. Format of refData. Accepts FASTA and other
               Biopython compatible formats. Optional if format is FASTA.
        """
        if name is None:
            print os.path.basename(self.resultsFile)
            name = "{}_PILEUP.txt".format(os.path.splitext(os.path.basename(self.resultsFile))[0])
        
        if sequences is None:
            sequences = self.records

        with open(name, "w") as save_file:
            for seq in sequences:
                new_seq = SeqRecord(Seq(str(seq.seq)), id=seq.id, description=seq.description)
                SeqIO.write(new_seq, save_file, file_format)

    def printTotalGaps(self, name=None):
        """Prints the global gap list to the standard out
        """
        if name is None:
           name = "Gaps_{}.txt".format(os.path.splitext(os.path.basename(self.resultsFile))[0])
        
        with open(name, "w") as gaps_file:
            print >> gaps_file, "index\tgap length"
            for index, gap in enumerate(self.total_gaps):
                print >> gaps_file, "{}\t{}".format(index, gap)

    def printCysCount(self, name=None):
        if name is None:
           name = "CysCount__{}.txt".format(os.path.splitext(os.path.basename(self.resultsFile))[0])
        
        with open(name, "w") as cys_file:
            print >> cys_file, "ranked_accession code\tnumber of Cys in match"
            for record in self.records:
                print >> cys_file, "{}\t{}".format(record.id, record.seq.count("C"))

class HMMSequence(MutableSeq):
    """An object that bridges an HMM hit sequence with Biopython's Seq object.
    This object determines the gaps of the hmm and insert them into the sequence
    and saves all of the information from the HMM hit into the sequence.

    Must subclass to work with different HMM profile packages."""
    
    def __init__(self, sequence, hmmLength, origSeqLength):
        """Intialise HMMSequence with the hmmer unit. Must run align and 
        determineGapPositions.

        Parameters:
        unit - HMMUnit object. 
        hmmLength - int. length of the HMM.
        align - boolean. Enables alignment algorithm based on HMM values.
                See HMMSequence.align() for more info.
        gaps - boolean. Enables gap counting algorithm to create pileup.
               See HMMSequence.determineGapPositions()
        """
        self.hmmLength = int(hmmLength)
        self.gaps = [0]*self.hmmLength
        self.origSeqLength = origSeqLength
        MutableSeq.__init__(self, sequence)

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
        for i in xrange(unknowns):
            self.data.insert(index, "X")

    def skip(self, percentX):
        """Returns true if the number X's is than than or equal to a given fraction 
        of the total hmm length and the number of consecuative X's is less than 
        to half of the given fraction of the hmm length

        Parameters:
        -----------
        percentX: float
            Percent Of total sequence is an X character
        """
        skip = self.count("X") > ceil(percentX*len(self)) and self.count("X"*int(floor(.5*percentX*self.hmmLength)))
        print "SKIP:", skip
        return skip

    def insertDeletionsAfter(self, index, deletions):
        """
        """
        for i, j in enumerate(range(deletions)):
            self.insert(index, "-")

    def printGaps(self, gaps=None):
        if gaps is None:
            gaps = self.gaps
        print "index\tgap length"
        for index, gap in enumerate(gaps):
            if gap == 0: continue
            print str(index) + "\t" + str(gap)

    def insertAllGaps(self, total_gaps):
        """This method inserts the previously computed gaps starting at the
        end of the sequence by reversing the list and recomputing.
        """
        total_lower = sum(self.gaps)
        refIndex = self.hmmLength-1
        index = len(self)-1
        while index >=0:
            refIndex = index-(total_lower-sum(self.gaps[refIndex:]))
            index -= self.gaps[refIndex]
            gap = total_gaps[refIndex]-self.gaps[refIndex]
            if gap>0:
                self.insertDeletionsAfter(index, gap)
            index-=1

    def align(self):
        """
        """
        raise NotImplementedError

    def determineGapPositions(self):
        """Determine how many gap positions and the indexes
        of the gaps in this hit. Must be called after 
        HMMSequence.aling() for accurate reading based on inseted X's
        """
        lowers = 0
        otherLetters = 0
        for i, c in enumerate(str(self)):
            if c.islower():
                lowers+=1
            else:
                otherLetters+=1
            if not c.islower() and lowers>0:
                gap_index = otherLetters-1
                self.gaps[gap_index] = lowers
                lowers = 0
