# HMMPileUP

Perform HMMER and UCSC-SAM model mixtures, saving the the results into a pile up. The sequences that have scores above a given thershold are dispayed with the portions that aligned to the HMM. Unknowns are filled in with X's in the front and end if sequence started or ended in the middle of the hmm. Next, gaps are calculated for all positions of the hmm in each sequence, and are reinserted into all sequences as dashes with appropriate lengths to make sure they are aligned. This work began in Dietlind L. Gerloff's (Ffame/UCSC) and is being continued for my PhD with Phil E. Bourne (NIH/NLM/NCBI) and Michael E. Grigg (NIH/UCSC).

## Requirements
* python2.7
* BioPython
* [HMMER 3.1.b2](http://hmmer.janelia.org/)
* [UCSC-SAM v3.5](https://compbio.soe.ucsc.edu/sam.html) (Note: convert.pl is required)

## Usage
1) Copy the Makefile into desired location and enter the locations of your reference data, target data, hmmer and sam binaries, and parameters for each model.

2) In your desired location, run:

`make target_dir`

where target is the name of the target data in fasta format without the file extension. The directory to find the target data must be specified in the Makefile.

3) Make a new directory and run HMMPileUP in that directory:

`python path/to/HMMPileUP/run.py reference_data --results results_directory`

Optionally, you can filter results to ignore and sequences that have too many X's by giving the "-X" or "--percentX" parameter with a percentage of  allowable X's compared to the length of the hmm.

## TODO
* Run pipeline from within run.py (Eliminates steps 1 and 2)
* Add examples

## Credits
* Dietlind L. Gerloff
* Edward Liaw (@edliaw) - Author of Makefile/pipeline
* Felicia Kemp
* Jonathan Magasin
