#! /usr/bin/env python
#PHE MLST software (PMS)  is a modified version of SRST version 1 script (http://sourceforge.net/projects/srst/files/?source=navbar), 
#modification made by  Anthony.Underwood@phe.gov.uk and Rediat.Tewolde@phe.gov.uk.
import os
import os.path
import sys
import argparse
import glob
import inspect
import pprint

module_folder_paths = ["modules"]
for module_folder_path in module_folder_paths:
    module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],module_folder_path)))
    if module_folder not in sys.path:
        sys.path.insert(1, module_folder)
from utility_functions import *
import log_writer
import MLST_extract_flanking_region_functions
import MLST_create_a_pileup_file
import MLST_read_pileup_file
import MLST_determine_score
import MLST_determine_ST
import MLST_report_ST_and_score

"""
This function checks if a file exists or not.

:param filepath: the name of the file
:param file_description: a description of the file that shown the error

:returns: error message that the file does not exist.
"""



def check_file_exists(filepath, file_description):
	if not os.path.exists(filepath):
		print("The " + file_description + " (" + filepath + ") does not exist")
		sys.exit(inspect.currentframe().f_lineno)
	else:
		print file_description + " detected: " + filepath

def parse_args(args):
	"""
	We have set parser = argparse.ArgumentParser() and added all arguments by adding parser.add_argument.
	"""
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-input_directory',
						'-i',
						help='please provide an input directory')
	parser.add_argument('--workflow',
						'-w',
						help='If using a workflow you must specify an input directory')
	parser.add_argument('--fastq_1',
						'-1',
						help='Fastq file pair 1')
	parser.add_argument('--fastq_2',
						'-2',
						help='Fastq file pair 2')
	parser.add_argument('--profile_file_directory',
						'-st',
						help='mlst directory')
	parser.add_argument('--output_directory',
						'-o',
						help='please provide an output directory')
	parser.add_argument('--bowtie',
						'-b',
						help='please provide the path for bowtie2', default='bowtie2')
	parser.add_argument('--samtools',
						'-sam',
						help='please provide the path for samtools', default='samtools')
	parser.add_argument('--log_directory',
						'-log',
						help='please provide the path for log directory')
	parser.add_argument('--infer_serotype',
						'-serotype',
						help='If you want to  infer salmonella serotype, please provide a True value')
	opts = parser.parse_args(args)
	return parser,opts


 
def main(parser, opts):
	
	fastq_files = []
	glob_pattern = "*fastq*"
	ids = None
	workflow_name = None
	version = None
	infer_salmonella_serotype = False
	
	# option 1: if user chooses to provide the workflow 
	if opts.workflow:
		
		#input_directory
		if not opts.input_directory:
			print("If using a workflow you must specify an input directory\n\n")
			print parser.print_help()
			sys.exit(inspect.currentframe().f_lineno)
		else:
			check_file_exists(opts.input_directory, 'input directory')
		
		#profile_file_directory
		profile_file_directory = os.path.dirname(os.path.realpath(__file__)) + "/MLST_data/" + opts.workflow
			
		#fastq_files
		fastq_files = glob.glob(opts.input_directory + "/" + glob_pattern)	
		if len(fastq_files) < 2:
			print("Fastq files are not pairs!")
			sys.exit(inspect.currentframe().f_lineno)
		
		#output_directory	
		if not opts.output_directory:
			opts.output_directory = opts.input_directory + '/mlst_typing'
			if not os.path.isdir(opts.output_directory): os.makedirs(opts.output_directory)
		
		#workflow_name
		workflow_name = opts.workflow
		
		#log_directory
		if not opts.log_directory:
			opts.log_directory = opts.input_directory + '/logs'
			if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory)
		
	##option 2: if user chooses to provide forward and reverse fastq_files then they can specify them with -1 and -2 options.
	elif opts.fastq_1 or opts.fastq_2:
		
		#fastq file
		check_file_exists(opts.fastq_1, 'Fastq 1')
		check_file_exists(opts.fastq_2, 'Fastq 2')
		fastq_files.append(opts.fastq_1)
		fastq_files.append(opts.fastq_2)
		

		
		#profile_file_directory
		if not opts.profile_file_directory:
			print("If you are passing full fastq paths to the fastq-1 and fastq-2 params, you need also to pass to the -st|--profile_file_directory  param a path to a reference_dir containing reference fasta files for each gene of interest and a profiles.txt file.txt')")
			print parser.print_help()
			sys.exit(inspect.currentframe().f_lineno)
		else:
			check_file_exists(opts.profile_file_directory, 'profile_file_directory')
			
		profile_file_directory = opts.profile_file_directory
		
		#output_directory
		if not opts.output_directory:
			output_directory = os.path.dirname(opts.fastq_1)
			opts.output_directory = output_directory + '/mlst_typing'
			if not os.path.isdir(opts.output_directory): os.makedirs(opts.output_directory)
		else:
			os.makedirs(opts.output_directory)
		
		#workflow_name
		workflow_name = "research"
		
		#log_directory
		if not opts.log_directory:
			opts.log_directory = opts.output_directory + '/logs'
			if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory) 
		
	
	if not os.path.exists(opts.output_directory + '/tmp'):
		os.makedirs(opts.output_directory + '/tmp')#make tmp directory in output_directory
	tmp_dir  = opts.output_directory + '/tmp'
	stderr_log_output = opts.log_directory + "/" + 'mlst_typing'+ ".stderr"
	stdout_log_output = opts.log_directory + "/" + 'mlst_typing'+ ".stdout"
	logger = log_writer.setup_logger(stdout_log_output, stderr_log_output)
	
	if opts.infer_serotype:
		infer_salmonella_serotype = True
		
	(SeqDir,seqFileName) = os.path.split(fastq_files[0])	
	(ids,ext) = os.path.splitext(seqFileName)
	version = "1-0"
	
	"""
	Extract flanking regions of 100bp upstream and downstream of each MLST locus by blast against a reference genome
	and concatenate flanking regions to correspondent locus variants sequence in fasta format.
	Newly concatenated sequence are then indexed by Bowtie2
	"""
	MLST_extract_flanking_region_functions.extract_flanking_region(tmp_dir,
															opts.output_directory,
															profile_file_directory,
															opts.bowtie,
															logger,
															stderr_log_output)
	

	"""
	Map each read set to each of the possible locus variants  by calling Bowtie2 (with  very sensitive options)
	and create tmp file. Convert the tmp to sam file by unset the secondary alignment bit score,
	convert the sam to BAM file, sort BAM file and then generate pileup file by using SAMtools mpileup.
	"""
	MLST_create_a_pileup_file.create_pileup_file(tmp_dir,
												 fastq_files,
												 opts.bowtie,
												 opts.samtools,
												 ids,
												 logger,
												 stderr_log_output)
	
	"""Parse through the pileup file to capture Depth of Coverage"""
	results, database, locusList, list_of_all_allele_numbers_tuple = MLST_read_pileup_file.get_scores(tmp_dir,
																									  profile_file_directory,
																									  stderr_log_output)
	
	"""Filter the correct allele based on  degree of variability of the read from the locus variant
	(by identifying present or absence of SNPs/INDELs)."""
	FinalResults= MLST_determine_score.determine_score(results,
													   locusList,
													   tmp_dir,
													   fastq_files,
													   opts.bowtie,
													   opts.samtools,
													   ids,
													   logger,
													   stderr_log_output)
	
	"""Determine ST values and write coverage  values and ST value to MLST_result.csv and results.xml files"""
	st = MLST_determine_ST.DetermineST(stderr_log_output,
									   database,
									   locusList,
									   list_of_all_allele_numbers_tuple,
									   FinalResults)
	
	MLST_report_ST_and_score.print_result(stderr_log_output,
										  st,
										  tmp_dir,
										  opts.output_directory,
										  profile_file_directory,
										  ids,
										  workflow_name,
										  version,
										  locusList,
										  FinalResults,
										  infer_salmonella_serotype)
	
	"""write ComponentComplete.txt file to the output directory"""
	try_and_except(stderr_log_output,
				   write_component_complete,
				   opts.output_directory)

   


if __name__ == "__main__":
    parser, opts= parse_args(sys.argv[1:])
    main(parser,opts)
