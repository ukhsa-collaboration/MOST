
import os
import sys
import pickle
import subprocess
import pprint
import Bio.SeqIO
from collections import Counter
import log_writer
from utility_functions import *
from MLST_create_a_pileup_file import create_bam_file	


def determine_score(results, locusList, tmp_dir, fastq_files, bowtie,samtools, ids, logger, stderr_log_output):	
	
	"""
	Function
	Filter the correct allele 

	The option for method:	
	LocusList[list]: list of locus names
	tmp_dir[str]: the path to  tmp directory
	fastq_files[str]: the path to the fastq file location
	bowtie[str]: the path to Bowtie2 command
	samtools[str]: the path to SAMtools command
	ids[str]: sample unique identifier number
	
	Return
	FinalResults[dict]
	"""

	FinalResults = try_and_except(stderr_log_output,
								  populate_scores_and_getNovelAllele,
								  results,
								  locusList,
								  tmp_dir,
								  fastq_files,
								  bowtie,
								  samtools,
								  ids,
								  logger)
	
	return FinalResults


def populate_scores_and_getNovelAllele(results, locusList, tmp_dir, fastq_files,  bowtie, samtools, ids, logger):
	
	"""
	Function
	Filter the correct allele
	
	The option for method:
	LocusList[list]: list of locus names
	tmp_dir[str]: the path to  tmp directory
	fastq_files[str]: the path to the fastq file location
	bowtie[str]: the path to  Bowtie2 command
	samtools[str]: the path to SAMtools command
	ids[str]: sample unique identifier number

	Return
	FinalResults[dict]
	"""
	
	ranges = pickle.load(open(os.path.join(tmp_dir, "ranges.pkl")))
	FinalResults = {}
	
	for locus in locusList:
		FinalResults[locus] = {}
		if locus not in results:
			FinalResults[locus]['VariantNumberHash'] ='-'
			FinalResults[locus]['ClosestVariantNumber']= '-'
			FinalResults[locus]['ReportedVariantNumber']='0'
			FinalResults[locus]['numberOfSNPs']= 0
			FinalResults[locus]['SNPsListsHash']= []
			FinalResults[locus]['numberOfINDELs']= 0
			FinalResults[locus]['INDELsListsHash']= []
			FinalResults[locus]['CoverageStat']= ("0/0/0/0/0/0/0")
			FinalResults = assign_CoverageStat(FinalResults,locus)
			FinalResults[locus]['percentage_coverage'] = 0
			continue
		
		vec = results[locus]  
		ExactMatch = filter(lambda itm: itm[2] == 0 and itm[7] == 0 and itm[10] == 100 , vec) 
		ExactMatch =sorted(ExactMatch, key=lambda tup:tup[8], reverse= True)
		
		#One significant matche with zero -SNP matches 
		if len(ExactMatch) > 0:
			FinalResults[locus]['VariantNumberHash'] = str(ExactMatch[0][1])
			FinalResults[locus]['ClosestVariantNumber']= str(ExactMatch[0][1]) 
			FinalResults[locus]['ReportedVariantNumber'] =str(ExactMatch[0][1])
			FinalResults[locus]['numberOfSNPs']= ExactMatch[0][2]
			FinalResults[locus]['SNPsListsHash']= ExactMatch[0][4]
			FinalResults[locus]['numberOfINDELs']= ExactMatch[0][7]
			FinalResults[locus]['INDELsListsHash']= ExactMatch[0][6]
			FinalResults[locus]['CoverageStat']= ('/'.join([str(i) for i in ExactMatch[0][5]]))
			FinalResults = assign_CoverageStat(FinalResults,locus)
			LocusNameWithVariantNumber = str(ExactMatch[0][3])
			LocusSize =  int(ranges[LocusNameWithVariantNumber][1]) - int(ranges[LocusNameWithVariantNumber][0])
			FinalResults[locus]['percentage_coverage'] = int(ExactMatch[0][8]/float(LocusSize)*100)

		#Multiple significant matches with zero -SNP matches 	
		if len(ExactMatch) > 1: 
			print  "locus " + ExactMatch[0][0] + "has multiple significant matches:"
			for e in ExactMatch:
				print  "\t" + e[1] 
		
		# No significant match, locus variant cannot be assigned. The closest variant is recorded  and the number of mismatches reported.	
		if len(ExactMatch) == 0:
			check_if_alleles_match_100_percent = filter(lambda itm: itm[10] == 100 , vec) # filter alleles which match 100% coverage
			if len(check_if_alleles_match_100_percent) > 0: # if allele match 100%, sort the filtered data by number of SNPs and extract the smallest allele varient
				vec = filter(lambda itm: itm[10] == 100 , vec)
				vec =sorted(vec, key=lambda tup:(tup[2], int(tup[1]))) # sort based on number of SNPs and allele varient
			else:
				vec =sorted(vec, key=lambda tup:(tup[2], int(tup[1]))) # sort based on number of SNPs and allele varient
			
			FinalResults[locus]['VariantNumberHash'] ='-'
			FinalResults[locus]['ClosestVariantNumber']= str(vec[0][1])
			FinalResults[locus]['ReportedVariantNumber'] ='*'+str(vec[0][1])
			FinalResults[locus]['numberOfSNPs']= vec[0][2]
			FinalResults[locus]['SNPsListsHash']= vec[0][4]
			FinalResults[locus]['numberOfINDELs']= vec[0][7]
			FinalResults[locus]['INDELsListsHash']= vec[0][6]
			FinalResults[locus]['CoverageStat']= ('/'.join([str(i) for i in vec[0][5]]))
			FinalResults = assign_CoverageStat(FinalResults,locus)
			LocusNameWithVariantNumber= str(vec[0][3])
			LocusSize =  int(ranges[LocusNameWithVariantNumber][1]) - int(ranges[LocusNameWithVariantNumber][0])
			FinalResults[locus]['percentage_coverage'] = int(vec[0][8]/float(LocusSize)*100)
			getNovelAllele(str(vec[0][1]),locus, fastq_files,  bowtie, samtools,ids,tmp_dir,logger) 
			
	return FinalResults


def assign_CoverageStat(FinalResults,locus):
	
	EachCoverageStat = FinalResults[locus]['CoverageStat']
	EachCoverageStat = EachCoverageStat.split('/')
	FinalResults[locus]['max_percentage_of_non_consensus_bases'] =EachCoverageStat[0]
	FinalResults[locus]['minimum_total_depth'] =EachCoverageStat[1]
	FinalResults[locus]['maximum_total_depth'] =EachCoverageStat[2]
	FinalResults[locus]['minimum_consensus_depth'] =EachCoverageStat[3]
	FinalResults[locus]['maximum_consensus_depth'] =EachCoverageStat[4]
	FinalResults[locus]['mean_consensus_depth'] =EachCoverageStat[5]
	
	return FinalResults
	

def getNovelAllele(variant, locus, fastq_files,  bowtie, samtools, ids, tmp_dir, logger):
	
	"""
	Function 
	Generate SAM, BAM and pileup file for novel allele
	
	The option for method:
	variant[str]: locus variant number
	locus[str]: locus name
	fastq_files[str]: the path to the fastq file 
	bowtie[str]: the path to  Bowtie2 command 
	samtools[str]: the path to SAMtools command
	ids[str]: sample unique identifier number
	tmp_dir[str]: the path to where  SAM, BAM, Pileup will be created
	"""
	
	refSeqs = pickle.load(open(os.path.join(tmp_dir, "refSeqs.pkl")))
	allele_name = locus + "-" + variant 
	typeFn = os.path.join(tmp_dir, allele_name + ".fa")
	typeFile = open(typeFn, "w")
	s = refSeqs[allele_name]
	Bio.SeqIO.write([Bio.SeqRecord.SeqRecord(s, id=allele_name)], typeFile, "fasta")
	typeFile.close()
	
	
	#index refrence sample
	bowtie2_index = bowtie + "-build"
	log_writer.info_header(logger, "index refrence sample")
	process = subprocess.Popen([bowtie2_index,
								typeFn,
								typeFn],
							stderr=subprocess.PIPE, stdout=subprocess.PIPE)## index ref
	process.wait()	
	log_writer.log_process(logger, process, log_error_to = "info")
	
	# create sam and bam files
	log_writer.info_header(logger, "Creating sam and bam files")
	bam = create_bam_file(tmp_dir,
						  fastq_files,
						  typeFn,
						  False,
						  bowtie,
						  samtools,
						  ids,
						  logger)
	
	#name bam file
	bamFn= os.path.join(tmp_dir, ids + "." + allele_name + ".bam")
	process = subprocess.Popen(['mv',bam,bamFn],stderr=subprocess.PIPE, stdout=subprocess.PIPE) 
	process.wait()
	log_writer.log_process(logger, process, log_error_to = "info")

	#index bam file
	process = subprocess.Popen([samtools, 'index',
								bamFn],
								stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	process.wait()
	log_writer.log_process(logger, process, log_error_to = "info")
	log_writer.info_header(logger, "index bam file")
	
	#generate pilup file
	piFn= os.path.join(tmp_dir,  ids + "." + allele_name + '.pileup')
	f = open(piFn, 'w')
	log_writer.info_header(logger, "generate pileup file")
	process = subprocess.Popen([samtools, 'mpileup',
								'-B', '-A', '-cf', typeFn, bamFn],
								stderr=subprocess.PIPE, stdout=subprocess.PIPE)# generate pileup
	for l in process.stdout:
		f.write(l)
	f.close()
	process.wait()	


