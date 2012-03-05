#!/usr/bin/python

from HMMERPileUp import HMMERPileUp
from SAMPileUp import SAMPileUp

import os, sys
import getopt

ref = ""
percentX = 0.5

def runHMMERPileUp(f):
    parser = HMMERPileUp(f)
    parser.parse(percentX=percentX)
    parser.addGapsToHMMSeqs()
    if not ref == "":
    	base = os.path.basename(f)
    	name = os.path.splitext(base)[0]
    	for refDB in os.listdir(ref):
    		if not refDB.endswith(".fasta"): continue
    		if os.path.splitext(refDB)[0] in name:
    			parser.setRefDB(os.path.join(ref, refDB))
    			break
    	if parser.refDB is None:
    		raise RuntimeError("No reference data at specified location")
    	parser.addGapsToRefData()
    parser.Write2File()
    
def runSAMPileUp(f):
    sam = SAMPileUp(os.path.splitext(f)[0], percentX=percentX)
    sam.Write2File()

def openDirectory(path, n=0):
	dir = os.listdir(path)
	for f in dir:
		print f
		if os.path.isdir(os.path.join(path,f)) and n<=2:
			openDirectory(os.path.join(path,f))
		if f.endswith("hmmer_hmmer.out") or f.endswith("sam_hmmer.out"):
			runHMMERPileUp(os.path.join(path,f))
			hmmer = True
		elif f.endswith("hmmer_sam.dist") or f.endswith("sam_sam.dist"):
			runSAMPileUp(os.path.join(path,f))
			sam = True

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["path=", "ref=", "percentX="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    path = os.getcwd()
    for opt, arg in opts:
        if opt in ["--path"]:
            if os.path.exists(arg):
                path = arg
            else:
            	raise RuntimeError("specified path doesn't exist")
        elif opt in ["--ref"]:
        	ref = arg
        elif opt in ["--percentX"]:
        	percentX = arg
        else:
            assert False, "unhandled option"
            
    openDirectory(path)


