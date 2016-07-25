# Performs HMMER/SAM model mixtures
# Edward Liaw 5/7/11

import os, sys
from subprocess import Popen

# Options
options = {
	"hmmer": {
		"build": {"--informat": "afa"},
		"search": {"--incdomE": 0.1}
	},
	"sam": {
		"build": {},
		"search": {"-select_mdalign":4, "-select_seq": 4, "-select_align": 4, "-select_score": 8, "-sort": 4, "-Emax": 10, "-mdEmax": 1}
	}
}

def run(reference_data, search_data, output=None, path=None):
	#Path to results files, keyesd by (train program, search program)
	hmmer = HMM("hmmer", reference_data, output=output, path=path)
	sam   = HMM("sam", reference_data, output=output, path=path)
	sam_hmmer = sam.convert()
	hmmer_sam = hmmer.convert()

	yield ("hmmer", "hmmer"), hmmer.search(reference_data)
	yield ("sam",   "hmmer"), sam_hmmer.search(reference_data)
	yield ("sam",   "sam"), sam.search(reference_data)
	yield ("hmmer", "sam"), hmmer_sam.search(reference_data)

class HMM(object):
	def __init__(self, program, reference_data, model_file=None, name=None, calibrate_sam=True, output=None, path=""):
		"""Train an HMM"""
		assert program.endswith("hmmer") or program.endswith("sam")
		self.program = program
		self.reference_data = reference_data
		self.model_file = model_file
		self.results_file = ""
		self.path = os.path.join(os.sep, "web", "public", "data", "projects", "6CysDB", "bin")
		self.output = output or os.getcwd()

		self.name = name or os.path.splitext(os.path.basename(reference_data))[0]

		#Train the HMM if no model file is given
		if self.model_file is None:
			self.model_file = os.path.join(
				self.output,
				"{}_{}".format(self.name, self.program)
			)
			if program.endswith("hmmer"):
				self.model_file += ".hmm"
				cmd = [os.path.join(self.path, "hmmbuild")]
				cmd += [str(item) for opt in options["hmmer"]["build"].items() for item in opt]
				cmd += [self.model_file, self.reference_data]
			elif program.endswith("sam"):
				self.model_file += ".mod"
				cmd = [os.path.join(self.path, "modelfromalign"), self.model_file[:-4], "-alignfile", self.reference_data]
				cmd += [str(item) for opt in options["sam"]["build"].items() for item in opt]
			else:
				raise RuntimeErrot("Invalid hmm name")

			print " ".join(cmd)
			process = Popen(cmd)
			process.communicate()
			process.wait()

		if self.program.endswith("sam") and calibrate_sam:
			self.calibrate()

	def search(self, search_data):
		""""""
		if self.name.endswith("hmmer"):
			#$(HMMSCORE) $(HMMOPTS) -o $@.out $^ $(APIDB)
			self.results_file = "{}.out".format(self.model_file[:-4])
			cmd = [os.path.join(self.path, "hmmsearch")]
			cmd += [str(item) for opt in options["hmmer"]["search"].items() for item in opt]
			cmd += ["-o", self.results_file, self.model_file, search_data]
			
		elif self.name.endswith("sam"):
			#$(SAMSCORE) $@ -i $^ -db $(APIDB) $(SAMOPTS)
			self.results_file = ["{}.{}".format(self.model_file, ext) \
				for ext in ("dist", "mstat", "mult")]
			cmd = [os.path.join(self.path, "hmmscore"), self.model_file[:-4], "-i", self.model_file, "-db", search_data]
			cmd += [str(item) for opt in options["sam"]["search"].items() for item in opt]

		process = Popen(cmd)
		process.communicate()
		process.wait()

		return self.results_file

	def convert(self):
		if self.program.endswith("hmmer"):
			#Convert hmmer to sam
			new_name = "hmmer2{}".format(self.program)
			new_model = "{}_{}.mod".format(self.reference_data, new_name)
			hmmer2_file = "{}.hmm".format(new_name)
			with open(hmmer2_file) as hmmer2:
				process = Popen([os.path.join(self.path, "hmmconvert"), "-2", self.model_file], stdout=hmmer2)
				process.wait()
			process = Popen([os.path.join(self.path, "convert"), hmmer2_file])
			process.wait()
			os.remove(hmmer2_file)
			#SAM convert renames file to have .con.asc.mod extension
			os.rename("{}.con.asc.mod".format(new_name), new_model)
		elif self.program.endswith("sam"):
			#Convert sam to hmmer
			new_name = "sam2{}".format(self.program)
			new_model = "{}_{}.hmm".format(self.reference_data, new_name)
			# Conversion script requires .asc.mod extension
			from shutil import copy
			sam_rename = "{}_{}.asc.mod".format(self.reference_data, new_name)
			copy(self.model_file, sam_rename)
			process = Popen([os.path.join(self.path, "hmmconvertSAM"), sam_rename])
			process.wait()
			os.remove(sam_rename)
			os.rename("{}.con.asc.hmm".format(new_name), new_model)
		else:
			raise RuntimeError("Must be a hmmer or sam model")

		return HMM(new_name, self.reference_data, model_file=new_model, path=self.path)

	def calibrate(self, num_sequences=1):
		"""Calibrate SAM models. HMMER automatically does this. From the SAM
		documentation"

		As can be seen, model calibration has added additional information to 
		the model library. The most important is the setting of the lambda 
		parameter for E-value calculation.

		If any of the parameters shown in the resulting model library are 
		changed, callibration must be performed again. It is important to note 
		that models must be calibrated separately for all differences in 
		scoring method.

		The process of model calibration involves scoring thousands of 
		sequences from either a database or a an internal random sequence 
		generator (for protein sequences, a dirichlet mixture distribution is 
		used). This can be a time-consuming process, but once calculated 
		hmmscore will produce much more accurate E-values.

		At present, we prefer scoring against a database of sequences. In this 
		process, we assume a symmetric distribution, and remove all reverse 
		null model scores less than (better than) zero to ensure that 
		calibration is based on non-matching sequences. To perform database 
		calibration, set calibrate to 1 to use the entire database, or a number 
		larger than one, such as 1000, to use the first 1000 sequences. The 
		evalues reported as a result of a calibration run will be the newly 
		calibrated evalues. If you generally score small sets of sequence, you 
		should calibrate your model against a large, non-redundent database, 
		and then score the calibrated model against the smaller database.

		For random sequence calibration, the only required parameter is 
		calibrate, which specifies the number of random sequences to score (if 
		calibrate is set to `1', an internal default is used). Optionally, 
		trackprior specifies a list (one per track) of Dirichlet mixtures over 
		sequence composition, genprot_prior indicates the default Dirichlet 
		mixture prior for protein sequences, genehl2_prior indicates the 
		default Dirichlet mixture prior for the EHL2 alphabet, gs_mean_log_len 
		is the natural logarithm of the mean sequence length to generate, and 
		gs_sd_log_len is the natural logarithm of the standard deviation of the 
		synthetic sequence length distribution. If a Dirichlet mixture is not 
		available for one or more of the sequence tracks, the characters are 
		drawn according to SAM's internal default background frequencies with 
		no variation in distribution between sequences.

		We are reasonably happy with our single-track model calibration method. 
		Our random sequence generators for multi-track model calibration are 
		not appropriately linked, and thus do not effectively calibrate 
		multi-track models, such as those with amino acid and secondary 
		structure sequences. This is currently an active research area, and we 
		hope to improve the calibration method for these models in the next 
		release.
		
		Parmeters
		---------
		num_sequences : int
			Number of random sequences to use from datbase. If 1, the entire database is used.
		"""
		process = Popen([os.path.join(self.path, "hmmscore"), self.name, "-modelfile", self.model_file, "-calibrate", str(num_sequences)])
		process.communicate()
