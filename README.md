#README
======
PHE MLST software (PMS) is a modified version of SRST (version 1) script.  
The purpose of this script is to assign MLST profiles and infer Salmonella serotyping to bacterial genomic sequence data.


##Table of content
----------------

  * Prerequisites
  * Instruction on how to download MLST database
  * Running PMS
  * Output


##Prerequisites
----------------
  * bowtie2/2.1.0
  * samtools/0.1.18 ( Please only use VERSION 0.1.18 )
  * emboss/6.6.0
  * blast+/2.2.27
  * python/2.7
  * yaml/1.1
  * numpy/python2.7/1.7.1
  * lxml/python2.7.0/3.2.3
  * biopython/python2.7/1.61


##Instruction on how to download MLST database
--------------------------------------------

1. Download profiles.txt file : profiles.txt file contains a list of STs and their allelic profiles.

2. Download loci variant sequences. One fasta file per locus saved as [name_of_locus].fas. The fasta sequence should be labelled with allele number
(e.g. aroc-1) and the locus label must be separated from the allele number by a dash ("-")
NB: There is a need to make sure the locus name exactly matches the labels used in the profiles.txt file.

3. Download a reference file in fasta format and save the file as reference.seq.  Subsequently, index
the reference sequence by using makeblastdb command. makeblastdb -in reference.seq -dbtype nucl  -out reference


##Usage
-------
	
	usage: PMS.py [-h]
	version 1-0, date 19/11/2015. 

PMS is a modified version of SRST version 1 script (http://sourceforge.net/projects/srst/files/?source=navbar, modification made by  Anthony.Underwood@phe.gov.uk and Rediat.Tewolde@phe.gov.uk.
	
	optional arguments:

   	-h, --help            show this help message and exit
  	-input_directory INPUT_DIRECTORY, -i INPUT_DIRECTORY
                        please provide an input directory
 	 --workflow WORKFLOW, -w WORKFLOW
                        If using a workflow you must specify an input
                        directory

  	--fastq_1 FASTQ_1, -1 FASTQ_1
                        Fastq file pair 1

  	--fastq_2 FASTQ_2, -2 FASTQ_2
                        Fastq file pair 2

  	--profile_file_directory PROFILE_FILE_DIRECTORY, -st PROFILE_FILE_DIRECTORY
                        MLST database

  	--output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
                        please provide an output directory

 	--bowtie BOWTIE, -b BOWTIE
                        please provide the path for bowtie2

 	--samtools SAMTOOLS, -sam SAMTOOLS
                        please provide the path for samtools

  	--log_directory LOG_DIRECTORY, -log LOG_DIRECTORY
                        please provide the path for log directory

  	--infer_serotype INFER_SEROTYPE, -serotype INFER_SEROTYPE
                        For Salmonella samples if you want to infer serotype,
                        please provide a True value
	
##Running PMS.py:

	A. To determine ST
	------------------------------
		Run the script by using the following 2 pathways:
		1.  The script can be run by specifying a workflow and full path to the input directory. The input directory should contain a pair of fastq files.
			e.g PMS.py -w <workflow name>  -i < full path to the input directory>
			Note:If you are using a workflow you must specify an input directory

		2. The script can be run by specifying the full path of the paired fastq files, MLST database and output directory
			e.g PMS.py -1 <full path to the fastq file> -2 <full path to the fastq file>  -o   <full path to the output directory>  -st <full path to the directory containing the alleles fasta and profiles.txt files>
			Note: If you are passing full fastq paths to the fastq -1 and fastq -2 params, you need also to pass to the -st|--profile_file_directory  param a path to a reference_dir containing 
		reference fasta files for each gene of interest and a profiles.txt files.


	B. To determine ST and infer salmonella serotype
	--------------------------------------------------
	Run the script by using the following 2 pathways:
	1.  The script can be run by specifying a workflow, full path to the input directory and set a serotype flag to True. The input directory should contain a pair of fastq files.
		e.g PMS.py -w <workflow name>   -i <path to the input directory> -serotype True
		Note:If using a workflow you must specify an input directory
		
	2. The script can be run by specifying the full path of the paired fastq files, MLST database, output directory and set a serotype flag to True.
		e.g PMS.py -1 <full path to the fastq file> -2 <full path to the fastq file> -o  <full path to the output directory>  -st <full path to the directory containing the alleles fasta and profiles.txt files> -serotype True
		Note: If you are passing full fastq paths to the fastq -1 and fastq -2 params, you need also to pass the -st|--profile_file_directory  param, a path to a reference_dir containing 
		reference fasta files for each gene of interest and a profiles.txt file.

 

##Output files
------------

1. <Sample ID>_MLST results.csv: describes allele designation and associated confidence quality metrics for each assignment

2. <sample ID>_results.xml: The XML file will report
	-	Sample ID : sample identifier
	-	Workflow value :  tested isolate, e.g. salmonella-typing
	-	Version : software version number 
	-	MLST value : ST value 
	-	Profile: allele variant for each locus
	-	QC mean consensus depth: Single value representing the minimum average consensus depth of all loci
	-	QC mean consensus depth of each locus: An average of the consensus depth across the full length of each locus. List of multiple values, one for each locus.
	-	QC max percentage non consensus base: The largest of the maximum percentage non consensus base values from all loci
	-	QC max percentage non consensus base of each locus: The maximum percentage non consensus base values across the full length of each locus. List of multiple values, one for each locus.
	-	QC minimum consensus depth: Single value representing the minimum consensus depth of all loci
	-	QC minimum consensus depth of each locus: Minimum consensus depth across the full length of each locus. List of multiple values, one for each locus.
	-	QC percentage coverage: Percentage coverage across allele length. NB: This value should always be 100% 
	-	Predicted_serotype: Predicted serotype and the number of occurrences associated with the ST in the PHE/Achtmann database. (for ONLY salmonella sample)
	-	Traffic light system: The script validates the results based on coverage metrics and writes a cut-off value standard based on the "Traffic light system". 
	The "Traffic light system" is assigned based on the following cut-off values:
		i.	The "GREEN traffic light" indication is assigned if the : 
			max percentage non consensus depth  < 15%  and
			Complete pileup= TRUE and
			Minimum consensus depth > 2 and
			Percentage coverage =100  and 
			ST not "Failed(incomplete locus coverage)

		ii.	The RED traffic light indication is assigned if the : 
			Complete pileup= FAIL or
 			Percentage coverage < 100 or
`			ST  is  "Failed(incomplete locus coverage)

		iii.	The "AMBER traffic light" indication is assigned if the there is no exact fit which matches either GREEN or RED


## Contact
-----------
Anthony.Underwood or Rediat Tewolde  
Anthony.Underwood@phe.gov.uk 
Rediat.Tewolde@phe.gov.uk
