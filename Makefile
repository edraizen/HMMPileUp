# Performs HMMER/SAM model mixtures
# Edward Liaw 5/7/11
# Folder paths
ROOT      = /Users/edraizen/Dropbox/Research/Gerloff/update2015/update_4_18_15
BIN       = /Users/edraizen/Dropbox/Research/Gerloff/bin
# Alignment file directory
DATA      = $(ROOT)

# Apicomplexa database path
APIDB     = $(ROOT)/EuPath.database2.4_18_15.fasta

# Program paths
HMMBUILD  = $(BIN)/hmmbuild
HMMSCORE  = $(BIN)/hmmsearch
SAMBUILD  = $(BIN)/modelfromalign
SAMSCORE  = $(BIN)/hmmscore
HMMCONV   = $(BIN)/hmmconvert
SAMCONV   = $(BIN)/convert
LOGO      = $(BIN)/makelogo
PDF       = ps2pdf12

# FASTA targets
TARGET   = PF04092_PF07422_seedmerged.Pf12D2_NMR_and_Pf12D1D2_Xray_MSAbasedUSETHIS.readyforapplication.141seqs_123pos.fasta \

# Options
FORMAT    = .fasta
HMMFORMAT = --informat afa
# For HMMER score
HMMOPTS   = --incdomE 0.1
# For SAM score
SAMOPTS   = -select_mdalign 4 -select_seq 4 -select_align 4 -select_score 8 -sort 4 -Emax 10 -mdEmax 1
# For logo maker
LOGOOPTS  = 

targets:
	for target in $(TARGETS); do \
		${MAKE} $$target\_dir ; \
	done

target:
	echo ${MAKE} $$(TARGET)\_dir
	${MAKE} $$(TARGET)\_dir

%_dir:
	echo "include ${CURDIR}/Makefile" > $*/Makefile	
	basename $* | xargs -0 -I name mkdir name
	cd $*;	${MAKE} $*_all
	
%_all:	
	${MAKE} $*_hmmer.hmm
	${MAKE} $*_sam.mod
	${MAKE} $*_hmmer_hmmer
	${MAKE} $*_hmmer_sam
	${MAKE} $*_sam_sam
	${MAKE} $*_sam_hmmer
	
%_hmmer_hmmer:	%_hmmer.hmm
	$(HMMSCORE) $(HMMOPTS) -o $@.out $^ $(APIDB)
	
%_hmmer_sam_nocal:	%_hmmer.hmm
	# Convert to HMMER 2.0 first
	$(HMMCONV) -2 $^ > $@.hmm
	$(SAMCONV) $@.hmm
	rm $@.hmm
	# Rename converted model
	mv $@.con.asc.mod $@.mod
	
	$(SAMSCORE) $@ -i $@.mod -db $(APIDB) $(SAMOPTS)
	
%_hmmer_sam:	%_hmmer.hmm
	# Convert to HMMER 2.0 first
	$(HMMCONV) -2 $^ > $@.hmm
	$(SAMCONV) $@.hmm
	rm $@.hmm
	# Rename converted model
	mv $@.con.asc.mod $@.mod
	${MAKE} $@.mlib
	
	$(SAMSCORE) $@ -i $@.mlib -db $(APIDB) $(SAMOPTS)
	
%_sam_hmmer:	%_sam.mod
	# Conversion script requires .asc.mod extension
	cp $^ $*_sam.asc.mod
	$(SAMCONV) $*_sam.asc.mod
	rm $*_sam.asc.mod
	# Rename converted model
	mv $*_sam.con.hmm $@.hmm
	
	$(HMMSCORE) $(HMMOPTS) -o $@.out $@.hmm $(APIDB)
	
%_sam_sam:	%_sam.mlib
	$(SAMSCORE) $@ -i $^ -db $(APIDB) $(SAMOPTS)
	
# Build HMMER model from FASTA file
%_hmmer.hmm:
	$(HMMBUILD) $(HMMFORMAT) $@ $(DATA)/$*$(FORMAT)
	${MAKE} $*_hmmer.pdf
	
# Build SAM model from FASTA file
%_sam.mod:
	$(SAMBUILD) $*_sam -alignfile $(DATA)/$*$(FORMAT)
	${MAKE} $*_sam.mlib
	${MAKE} $*_sam.pdf
	
# Calibrated SAM model
%.mlib:	%.mod
	$(SAMSCORE) $* -modelfile $^ -calibrate 1
	
# Make logos
%_hmmer.pdf:	%_hmmer.hmm
	$(HMMCONV) -2 $^ > $*_tmp.hmm
	$(SAMCONV) $*_tmp.hmm
	$(LOGO) $*_hmmer -modelfile $*_tmp.con.asc.mod $(LOGOOPTS)
	$(PDF) $*_hmmer.eps
	rm *_tmp*
	
%_sam.pdf:	%_sam.mod
	$(LOGO) $*_sam -modelfile $^ $(LOGOOPTS)
	$(PDF) $*_sam.eps
	
