from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
import sys
import re
            
class Sequences(list):
    def __init__(self, seqFile, format="fasta"):
        for seq in SeqIO.parse(seqFile, format):
            self.append(seq)

    def __getattribute__(self, name):
        try:
            return list.__getattribute__(self, name)
        except AttributeError:
            #Make it so the user can easily get access to the first item, useful when they know there is only element
            return self[0].__getattribute__(name)
            
    def numberOfSequences(self):
        return len(self)

    def numberOfOrganisms(self):
         if not hasattr(self, "species"):
            self.getOrganisms()
         return len(self.species)

    def getOrganisms(self):
        self.organisms = []
        for seq in self:
            arrName = re.search("\[([^}]+)\]", seq.id).group()[1:-1]
            if "[" in arrName:
                arrName = arrName.split("[")[1]
            arrName = arrName.split(" ")
            name = "%s %s" % (arrName[0], arrName[1])
            if name is not None and name not in self.organisms and not "[" in name and not "]" in name:
                self.organisms.append(name)
                seq.organism = name
        self.species.sort()
        return self.species

    def setHeaders(self, fmt):
        if "{Sp}" in fmt and not hasattr(self, "species"):
            self.getOrganisms()
            
        for seq in self:
            header = fmt
            header.replace("{GI}", seq.id.split("|")[1])
            header.replace("{Sp}", seq.organism)
            header.replace("{Len}", len(seq))
            seq.id = header

    def sequenceLengthRange(self):
        seqs = self[:]
        seqs.sort(key=lambda x: len(x))
        print "HUGE: ", seqs[-1].comment
        self.lengths = [len(seq) for seq in seqs]
        return self.lengths[0], float(sum(self.lengths))/len(self.lengths), self.lengths[-1]

    def plotLengths(self):
        lengths = [len(seq) for seq in self.sequences]
        ind = np.arange(len(lengths))  # the x locations for the groups
        width = 0.35       # the width of the bars

        fig = plt.figure()
        ax = fig.add_subplot(111)
        rects1 = ax.bar(ind, lengths, width)
        ax.set_ylabel('Sequence Length')
        ax.set_xlabel('Apicomplexa FASTA refSeq sequences')
        ax.set_title('Sequence Length for Apicomplexa')
        plt.show() 
    
    def get(self, index=0):
        return self[index]
    
    def toFile(self, name):
        SeqIO.write(self, name, "fasta")

class MutableSequences(Sequences):
    def __init__(self, seqFile, format="fasta"):
        for seq in SeqIO.parse(seqFile, format):
            seq.seq = MutableSeq(seq.seq.tostring())
            self.append(seq)

class EntrezSequences(Sequences):
    def __init__(self, codes, email, db="protein", rettype="fasta"):
        print "code:", codes, codes.__class__.__name__
        Entrez.email = email
        try:
             f = open(codes)
             codes = []
             for code in f.readlines():
                 codes.append(code)
             f.close()
        except:
            pass
        print codes
        
        for code in codes:
            handle = Entrez.efetch(db=db, rettype=rettype, id=code)
            print handle.url
            seq_record = SeqIO.read(handle, "fasta")
            self.append(seq_record)
            handle.close()

if __name__ == "__main__":
    codes = EntrezSequences(sys.argv[1], "edraizen@ucsc.edu")
    print codes.numberOfOrganisms()
    ## for code in f:
    ##     codes.append(code);
    ## seqs = sequence(codes)
    ## seqs.parse()
    ## f.close()
    
