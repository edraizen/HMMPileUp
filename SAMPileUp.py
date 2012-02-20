import os, sys
import re
import getopt
from Bio.Seq import Seq, MutableSeq
from Bio import SeqIO
from Sequences import MutableSequences
from Bio import SeqIO
from time import sleep

class SAMPileUp(object):
    def __init__(self, inFile, evalcutoff=0.0):
        self.inFile = inFile
        self.distFile = "%s.dist" % inFile #Get HMM Length
        self.multFile = "%s.mult" % inFile #Get sequences
        self.hmmLength = None
        self.evalcutoff = evalcutoff

        #Import sequences from .mult file
        self.sequences = [] #MutableSequences(self.multFile)
        for seq in SeqIO.parse(self.multFile, "fasta"):
            self.sequences.append(seq)

        dist = open(self.distFile, "r")
        getEvalues = False
        for i, line in enumerate(dist.readlines()):
            try:
                if line.split(",")[0].split(" ")[-1] == "sequences":
                    self.hmmLength = int(line.split(", ")[2].split(" ")[0])
                else:
                    pass
            except:
                continue
            if "% Sequence ID" in line:
            	getEvalues = True
            	print "Assigning evaules to sequences"
                continue
            if getEvalues:
                sys.stdout.write("\r" + "."*(i%4))
                sys.stdout.flush()
                try:
                    _id, length, simple, reversed, evalue = line.split()[:5]
                except:
                    #raise RuntimeError("Invalid mult file")
                    break
                if float(evalue)<self.evalcutoff:
                    continue
                for seq in self.sequences:
                    seq.evalue = -2
                    if _id in seq.id:
                        seq.description = "[score: %s bits]" % evalue
                        seq.evalue = float(evalue)
                        break
        print ""               
                        
        dist.close()

        if self.hmmLength is None:
            assert "Error! %s is not valid" % self.distFile
                  
        self.align()
        self.insertAllGaps()

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
	
		self.gap=[0]*self.hmmLength
		removedSequences = []
		
		for idx, record in enumerate(self.sequences):
			seq = record.seq.tostring()
			m = re.search('([a-z]*)(-*)([A-Z][A-Za-z\-]*[A-Z])(-*)([a-z]*)', seq)
			frontLowers = m.groups(0)[0]
			frontDashes = m.groups(1)[1]
			matchedSeq = m.groups(1)[2]
			endDashes = m.groups(1)[3]
			endLowers = m.groups(1)[4]
			testSeq = "%s%s%s" % (frontDashes, matchedSeq, endDashes)
		  
			#Append X's to beginning
			if len(frontDashes)<=len(frontLowers):
			    realSeq = "X"*len(frontDashes)
			else:
				realSeq = "-"*(len(frontDashes)-len(frontLowers))
				realSeq += "X"*len(frontLowers)
					
			realSeq += matchedSeq
	
			#Append X's to end
			if len(endDashes)<=len(endLowers):
				realSeq += "X"*len(endDashes)
			else:
				realSeq += "X"*len(endLowers)
				realSeq += "-"*(len(endDashes)-len(endLowers))
	
			lowers = 0
			otherLetters = 0
			record.gap = [0]*self.hmmLength
			for i, c in enumerate(realSeq):
				if c.islower():
					lowers+=1
				else:
					otherLetters+=1
				if not c.islower() and lowers>0:
					gap_index = (otherLetters-1)
					record.gap[gap_index] = lowers
					if self.gap[gap_index] < lowers:
						self.gap[gap_index] = lowers
					lowers = 0      
					
			#Remove sequence if there is more 50% X in the string
			if record.seq.count("X")>len(record.seq)/2:
				removedSequences.append(record.decription)
				del self.sequences[idx]
				
			record.seq = MutableSeq(realSeq)
			record.seq.count
		if len(removedSequences) > 0:
			print "Removed sequences from %s:" % os.path.basename(self.inFile)
			for seq in removedSequences:
				print seq
		#self.printGaps()

    def printGaps(self):
        print "index\tgap length"
        for index, gap in enumerate(self.gap):
            if gap == 0: continue
            print str(index) + "\t" + str(gap)

    def insertAllGaps(self):
        """This method inserts the previously computed gaps starting at the
        end of the sequence by reversing the list and recomputing.
        """
        for record in self.sequences:
            total_lower = sum(record.gap)
            refIndex = self.hmmLength-1
            index = len(record.seq.data)-1
            while index >=0:
                refIndex = index-(total_lower-sum(record.gap[refIndex:]))
                index -= record.gap[refIndex]
                gap = self.gap[refIndex]-record.gap[refIndex]
                if gap>0:
                    for i, j in enumerate(range(gap)):
                        record.seq.data.insert(index, "-")
                index-=1
                        
    def Write2File(self, name=None, type="fasta"):
        """Save the hit sequences to any format
        name - string. Name to save the ref data, if not same name.
        type - string. Format of refData. Accepts FASTA and other
        Biopython compatible formats. Optional if format is FASTA.
        """
        if name is None:
            name = "%s_PILEUP.txt" % (os.path.basename(self.inFile))
        
        if name == "PF07422_seed_PILEUP.txt" or name == "PF04092_seed_PILEUP.txt" or os.path.exists(name):
        	print "File name (%s) already exists. %s, %s"%(os.path.basename(self.inFile), os.path.splitext(os.path.basename(self.inFile))[0], name) 
        	name=raw_input("Pick a new name:")
        
        tmp = []
        for i, seq in enumerate(self.sequences):
            tmp_seq = seq
            if not seq.__class__.__name__ == "SeqRecord":
                tmp_seq.seq = seq.seq.toseq()
            tmp_seq.id = "%d_%s"%(i+1, tmp_seq.id)
            tmp.append(tmp_seq)
        SeqIO.write(tmp, name, type)

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["eVal="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    evaluecutoff = 0.0
    for opt, arg in opts:
        if opt in ["--eVal"]:
            evaluecutoff = float(arg)
        else:
            assert False, "unhandled option"
            
    sam = SAMPileUp(sys.argv[1], evaluecutoff)
    sam.Write2File()            

            

