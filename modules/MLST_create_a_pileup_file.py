
import os
import os.path
import sys
import subprocess
import fileinput
import log_writer
from utility_functions import *


def create_pileup_file(tmp_dir, fastq_files, bowtie, samtools, ids, logger, stderr_log_output):
	
	"""
	Function
	(1) Map each read set to each of the possible locus variants  by calling Bowtie2
	(with  very sensitive options) and create  tmp file
	(2) Convert the tmp to sam file by unset the secondary alignment bit score
	(3) Convert the sam to BAM file
	(4) Sort BAM file
	(5) Generate pileup file by using SAMtools mpileup command
	
	The option for method:
	tmp_dir[str]: the path to where the tmp, SAM, BAM and sorted BAM files will be created
	fastq_files[list]: the path to the paired fastq files location
	bowtie[str]: the path to Bowtie2 command
	samtools[str]: the path to SAMtools command
	ids[str]: sample unique identifier number
	"""
	
	refFn = os.path.join(tmp_dir, "reference.fa")
	expand = True
	sorted_bam_file = try_and_except(stderr_log_output,
									 create_bam_file,
									 tmp_dir,
									 fastq_files,
									 refFn,
									 expand,
									 bowtie,
									 samtools,
									 ids,
									 logger)
	try_and_except(stderr_log_output,
				   pileupReads,
				   tmp_dir,
				   sorted_bam_file,
				   refFn,
				   samtools,
				   logger)


def create_bam_file(tmp_dir, fastq_files, refFn, expand, bowtie, samtools, ids, logger):
	
	"""
	Function
	(1) Map each read set to each of the possible locus variants  by calling Bowtie2
	(with  very sensitive options) and create tmp file
	(2) Convert the tmp to sam file by unset the secondary alignment bit score
	(3) Convert the sam to BAM file
	(4) sort BAM file
	
	The option for method:
	tmp_dir[str]: the path to where the tmp, SAM, BAM and sorted BAM files will be created
	fastq_files[list]: the path to the fastq file location
	refFn[str]:  the path to the  reference file location
	expand[str] : True
	bowtie[str]: the path to  Bowtie2 command 
	samtools[str]: the path to SAMtools command
	ids[str]: unique identifier number
	logger[str]: the path to where the stderr and stdout logged
	
	Return
	sorted_bam_file[str]: sorted BAM file
	"""
	
	tmp = os.path.join(tmp_dir, ids +'.tmp') # temporary sam output
	sam = os.path.join(tmp_dir, ids +  '.sam')
	
	if expand:
		#1. Creating tmp file
		log_writer.info_header(logger, "Creating tmp file")
		# -k = report up to 99999 good alignments per read.
		#--very-sensitive option = -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 
		process = subprocess.Popen([bowtie, '--fr',
									'--no-unal',
									'--minins', '300', '--maxins', '1100',
									'-x', refFn,
									'-1', fastq_files[0], '-2', fastq_files[1],
									'-S', tmp,
									'-k', '99999', '-D', '20', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50'],
									stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		process.wait()
		log_writer.log_process(logger, process, log_error_to = "info")
		
		#2.remove_secondary_mapping_bit
		log_writer.info_header(logger, "remove_secondary_mapping_bit")
		i = open(tmp)
		o = open(sam, 'w')
		remove_secondary_mapping_bit(tmp, sam)
		i.close()
		o.close()
	
	else:
		log_writer.info_header(logger, "Creating sam file")
		process= subprocess.Popen([bowtie,  '--fr',
								   '--no-unal',
								   '--minins', '300', '--maxins', '1100',
								   '-x', refFn,
								   '-1', fastq_files[0], '-2', fastq_files[1],
								   '-S', sam,
								   '-k', '99999', '-D', '20', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50'],
			stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		process.wait()
		log_writer.log_process(logger, process, log_error_to = "info")
	
	#3.Converting sam to bam file
	bam = os.path.join(tmp_dir, ids +  '.unsortedbam')
	log_writer.info_header(logger, "Converting sam to bam")
	process = subprocess.Popen([samtools, 'view',
								'-bhS',
								'-o', bam, sam], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	process.wait()
	log_writer.log_process(logger, process, log_error_to = "info")
	
	#4.Sort bam file
	out0 = os.path.join(tmp_dir,ids +  '-all')
	log_writer.info_header(logger, "Sorting bam")
	sorted_bam_file = os.path.join(tmp_dir, ids + '-all.bam')
	process = subprocess.Popen([samtools, 'sort',
								bam, out0], stderr=subprocess.PIPE, stdout=subprocess.PIPE) 
	process.wait()
	log_writer.log_process(logger, process, log_error_to = "info")
	
	return sorted_bam_file


def remove_secondary_mapping_bit(sam, sam_parsed):
	
	"""
	Function
	Takes a SAM file and deducts 256 from the second column(FLAG) that unset the
	secondary alignment bit score 
	NB: reads with bit(250) set are not reported when using Samtools pileup
	
	The option for method:
	sam[string]: SAM file 
	sam_parsed[string]: parsed SAM file
	"""
	
	lines = iter(fileinput.input([sam]))
	sam_parsed_file = open(sam_parsed, "w")
	headers = []
	body = []
	
	for line in lines:
		if line.startswith('@'):
			sam_parsed_file.write(line)
		else:
			# chomp line
			line = line.rstrip('\n')
			details = line.split("\t")
			flag = int(details[1])
			if flag > 256:
				details[1] = str(flag - 256)
			print >> sam_parsed_file, '\t'.join(details)
	sam_parsed_file.close()


def pileupReads(tmp_dir, sorted_bam_file, refFn, samtools, logger):
	
	"""
	Function
	Generate pileup file by using SAMtools mpileup command.
	NB: use -B -A -f option to optimises coverage and --A
	flag count anomalous read 
	 
	The option for method:
	tmp_dir[str]: the path to where  pileup file will be created
	sorted_bam_file[str]:  the path to the BAM file location
	refFn[str]:  the path to the reference file location
	samtools[str]: the path to SAMtools command
	logger[str]: the path to where the stderr and stdout logged
	"""
	
	#1. Index bam file
	log_writer.info_header(logger, "index bam file")
	process = subprocess.Popen([samtools, 'index',
								sorted_bam_file],
								stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	process.wait()
	log_writer.log_process(logger, process, log_error_to = "info")
	
	#2. Generate pileup file
	pileFn =  os.path.join(tmp_dir, 'all.pileup')
	pileupFile = open(pileFn, 'w')
	log_writer.info_header(logger, "Generate pileup file")
	process = subprocess.Popen([samtools, 'mpileup', '-B', '-A', '-f',
								refFn, sorted_bam_file],
								stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	for l in process.stdout:
		pileupFile.write(l)
	process.wait()
	log_writer.log_process(logger, process, log_error_to = "info")
	pileupFile.close()
	
	
	