#Metric-Oriented Sequence Typer (MOST) software is a modified version of SRST version 1 script (http://sourceforge.net/projects/srst/files/?source=navbar), 
#modification made by  Anthony.Underwood@phe.gov.uk and Rediat.Tewolde@phe.gov.uk.
import  re
import os
import os.path
import sys
import subprocess
import pickle
import glob
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import Bio.Seq
import Bio.SeqIO
import log_writer
from utility_functions import *


def extract_flanking_region(tmp_dir, output_directory, profile_file_directory, bowtie, logger, stderr_log_output):
	
	"""
	Function
	(1) Extract flanking regions of 100bp upstream and downstream of each MLST locus by blast against a
	reference genome. BLAST uses the first locus sequence as a query.
	(2) Concatenate flanking regions to correspondent locus variants sequence in fasta format.
	Newly concatenated sequence are then indexed by Bowtie2
	(3) Then extract and store as pickled object:
		a. locus- variant names (loci.pkl)
		b. start and end position of locus variant sequences (without the flanking sequences)
		c. Locus variants sequence (refSeqs.pkl)
	
	The option of the method
	tmp_dir[str]: the path to where refSeqs.pkl, ranges.pkl and loci.pkl will be created
	output_directory[str]: the path to where the summary.txt file will be created
	profile_file_directory[str]: the path to the reference.seq, profile.txt and
	the  Locus variant sequences (*.fas) files location
	bowtie[str]: the command used to index the reference sequence
	"""
	
	try_and_except(stderr_log_output,
				   flanking_regions,
				   profile_file_directory,
				   output_directory,
				   logger)
	
	try_and_except(stderr_log_output,
				   concatenate_flanking_regions,
				   output_directory + "/summary.txt",
				   tmp_dir,
				   bowtie,
				   logger)

	
def flanking_regions(profile_file_directory, output_directory, logger):
	
	"""
	Function
	(1) Extract flanking regions of 100bp upstream and downstream of each MLST locus
	by blast against a reference genome.BLAST uses the first locus sequence as a query.
	(2) Creates summary.txt file (a tab-delimited text file display the path to the
	loci and flanking sequences) 
	
	The option of the method
	profile_file_directory[str]: The path to the reference.seq, profile.txt and
	the Locus variant sequences (*.fas) files location
	output_directory[str]: The path to where the summary.txt file will be created
	logger[str]: The path to where the stderr and stdout logged
	"""
	
	reference_fasta_file = profile_file_directory + "/reference.seq"	
	refseq_record = SeqIO.read(reference_fasta_file, "fasta", generic_dna)
	locus_files = glob.glob(profile_file_directory + "/*.fas")
	locus_files = sorted(locus_files)
	
	summary_file_handle = open(output_directory + "/summary.txt", "w")
	for seq in locus_files:
		(seqDir,seqFileName) = os.path.split(seq)	
		(seqBaseName,ext) = os.path.splitext(seqFileName)
		bait = seqBaseName + "_bait.fasta"
		log_writer.info_header(logger, "create bait file")
		process = subprocess.Popen(['seqret',seq,'-firstonly','-auto','-out',output_directory+ '/' + bait], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		process.wait()
		log_writer.log_process(logger, process, log_error_to = "info")
		cline = NcbiblastnCommandline(query=output_directory+ '/' + bait, db=profile_file_directory + "/reference",evalue=0.001, out=output_directory + "/my_blast_tmp.xml", outfmt=5)		
		stdout_log_output, stderr_log_output = cline()
		result_handle = open(output_directory + "/my_blast_tmp.xml")
		blast_record = NCBIXML.read(result_handle)
		query_length = blast_record.query_letters
		for alignment in blast_record.alignments:
			hsp = alignment.hsps[0] 
			if hsp.align_length/float(query_length) > 0.5:
				if hsp.sbjct_start > hsp.sbjct_end:
					subject_start = hsp.sbjct_start + (hsp.query_start - 1)
				else:
					subject_start = hsp.sbjct_start - (hsp.query_start - 1)
				if hsp.sbjct_start > hsp.sbjct_end:
					subject_end = hsp.sbjct_end - (query_length - hsp.query_end)
				else:
					subject_end = hsp.sbjct_end + (query_length - hsp.query_end)
				revcomp = 1 
				if hsp.sbjct_start > hsp.sbjct_end:
					revcomp = -1
				left_coords = [min(subject_start,subject_end)-100,min(subject_start,subject_end)-1]
				right_coords = [max(subject_start,subject_end)+1,max(subject_start,subject_end)+100]
				left_cmd = ["seqret ",reference_fasta_file," -sbegin ",str(left_coords[0])," -send ",str(left_coords[1])," -osformat fasta -auto -out " + output_directory + "/tmp_left_flank.fasta"]
				os.system(''.join(left_cmd))
				
				right_cmd = ["seqret ",reference_fasta_file," -sbegin ",str(right_coords[0])," -send ",str(right_coords[1])," -osformat fasta -auto -out " + output_directory + "/tmp_right_flank.fasta"]
				os.system(''.join(right_cmd)) 
				
				left_record = SeqIO.read(output_directory + "/tmp_left_flank.fasta", "fasta")
				
				if revcomp < 0:
					left_record.id = "down"
					left_record.seq = left_record.seq.reverse_complement()
				else:
					left_record.id = "up"
				right_record = SeqIO.read(output_directory + "/tmp_right_flank.fasta", "fasta")
				if revcomp < 0:
					right_record.id = "up"
					right_record.seq = right_record.seq.reverse_complement()
				else:
					right_record.id = "down"
				right_record.description = ""
				left_record.description = ""
				out_handle = open(output_directory + "/" + seqBaseName + "_flanks.fasta", "w")
				out_handle.write(right_record.format("fasta"))
				out_handle.write(left_record.format("fasta"))
				out_handle.close()
				summary_file_handle.write('\t'.join([seqBaseName,seq,output_directory + "/" + seqBaseName + "_flanks.fasta"]) + "\n")
	summary_file_handle.close()
	
				

def concatenate_flanking_regions(specFn, tmp_dir, bowtie, logger):
	
	"""
	Function
	(1) Concatenate flanking regions to correspondent locus variants sequence in fasta format.
	Newly concatenated sequence are then indexed by Bowtie2
	(2) Then extract and store as pickled object:
		a. locus- variant names (loci.pkl)
		b. start and end position of locus variant sequences (without the flanking sequences)(ranges.pkl)
		c. Locus variants sequence (refSeqs.pkl)
	
	The option of the method
	specFn[str]: A tab-delimited text file display the path to the seven flanking  and loci sequences(summary.txt)
	tmp_dir[str] The path to where refSeqs.pkl, ranges.pkl and loci.pkl will be created
	bowtie[str]: The command used to index the reference sequence
	logger[str]: The path to where the stderr and stdout logged
	
	Return
	loci[list]: loci name
	"""
	
	(specDir,summaryFileName) = os.path.split(specFn)
	spc = []
	
	for l in open(specFn):
		spc.append(l.split())
	refFn = os.path.join(tmp_dir, "reference.fa")
	rf = open(refFn, "w") 
	ranges = {}
	loci = [] 
	refSeqs = {} 
	for (loc, variantsFn, flanksFn) in spc:
		loci.append(loc)
		fs = {} 
		f = open(os.path.join(specDir, flanksFn))
		for r in Bio.SeqIO.parse(f, "fasta"):
			fs[r.id] = r.seq
		f = open(os.path.join(specDir, variantsFn))
		for r in Bio.SeqIO.parse(f, "fasta"):
			s = Bio.Seq.MutableSeq('', Bio.Alphabet.generic_dna) 
			s += fs['up']
			s += r.seq
			s += fs['down']
			Bio.SeqIO.write([Bio.SeqRecord.SeqRecord(s, id=r.id)], rf, "fasta") 
			ranges[r.id] = (len(fs['up']), len(fs['up']) + len(r.seq))
			refSeqs[r.id] = s 
	rf.close()
	rangesFn = os.path.join(tmp_dir, "ranges.pkl") #start and end position of locus variant sequences (without the flanking sequences)
	f = open(rangesFn, 'w')
	pickle.dump(ranges, f)
	f.close()
	lociFn = os.path.join(tmp_dir, "loci.pkl")
	f = open(lociFn, 'w')
	pickle.dump(loci, f)
	f.close()
	refSeqsFn = os.path.join(tmp_dir, "refSeqs.pkl") #Locus variants sequence
	f = open(refSeqsFn, 'w')
	pickle.dump(refSeqs, f)
	f.close()
	bowtie2_index = bowtie + "-build"
	log_writer.info_header(logger, "bowtie_indexed")
	process = subprocess.Popen([bowtie2_index, refFn, refFn], stderr=subprocess.PIPE, stdout=subprocess.PIPE) # generate index of reference fasta for mapping
	process.wait()
	
	log_writer.log_process(logger, process, log_error_to = "info")
	os.system("rm -f summary.txt")
	
	



