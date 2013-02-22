from Bio import SeqIO

class PileUpIOMixins(object):
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

    def printTotalGaps(self):
        """Prints the global gap list to the standard out
        """
        
        print "index\tgap length"
        for index, gap in enumerate(self.total_gaps):
            if gap == 0: continue
            print str(index) + "\t" + str(gap)
            
    def printCysCount(self):
    	out = "ranked_accession code\tnumber of Cys in match\n"
    	for record in self.records:
    		out += "%s\t%d\n" % (record.id, record.seq.count("C"))
    	f = open("CysCount_%s"%(os.path.splitext(os.path.basename(self.resultsFile))[0]), 'w')
    	f.write(out)
    	f.close()