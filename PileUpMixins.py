class PileUpMixins(object):
	def ReIndexSequences(self):
        for i, record in enumerate(self.records):
            record.id = "{}_{}".format(i+1, record.id.split("_", 1))

    def removeSeqs(self, percentX=0.5):
        #Remove sequence if there is more %percentX in the string
        #BUG: float conv didn't work, fixed, but commented out to test previous results
        minAAcount = percentX*self.hmmLength
        tmpRecords = []
        for record in self.records:
            print len(str(record.seq))
            nonSeqCharacters = str(record.seq).translate(None, "ACDEFGHIKLMNPQRSTVWY")
            if nonSeqCharacters < minAAcount:
                #Remove record, warn
                print "Removed sequences from %s:" % os.path.basename(self.resultsFile)
                print seq
            else:
                #Save record
                tmpRecords.append(record)

        self.records = tmpRecords
        self.ReIndexSequences()

    def ignoreBitsBelow(self, bitThreshhold):
        """Ignore bits below a certain threshhold.

        Parameters:
        bitThreshold - a floating point value. Ignore bits that are less than bitThreshold.
        """
        tmpRecords = []
        for record in self.records:
            if record.hmm.bits < threshold:
                tmpRecords.append(record)

        self.records = tmpRecords

    def addGapsToRefData(self, refDB, name=None, format="fasta"):
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
