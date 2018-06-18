#README
======
Metric-Oriented Sequence Typer (MOST) software is a modified version of SRST (version 1) script by Dr Kat Holt and colleagues (http://sourceforge.net/projects/srst/files/).
The purpose of this script is to assign MLST profiles and infer Salmonella serotyping from bacterial genomic short read sequence data.


##Table of content
----------------

  * Prerequisites
  * Instructions on how to download MLST database
  * Running MOST
  * Output


##Prerequisites
----------------
  * bowtie2 2.1.0 or later
  * samtools 0.1.18 ( Please only use VERSION 0.1.18 )
  * emboss 6.6.0
  * blast+ 2.2.27
  * python 2.7
  * yaml 1.1
  * numpy 1.7.1
  * lxml 3.2.3
  * biopython 1.61


##Instruction on how to download MLST database
--------------------------------------------

Download ST profiles and allele sequences from the MLST databases (e.g http://pubmlst.org/data/ , http://www.mlst.net/databases/ ) and reference complete genome

1. ST profiles: should be tab-delimited text file and should contain a list of STs and their allelic profiles. The file should be named as profiles.txt. Please note that profiles.txt header should follow the format given in MLST databases
(e.g. ST TAB followed by locus names). Any extra columns after the locus column should be removed ( e.g. clonal complex)

2. Loci variant sequences: should be in fasta format. The fasta file per locus should be named as [name_of_locus].fas (e.g. aroc.fas).
The fasta sequence headers should be labelled with their allele number(e.g. aroc-1) and the locus name must be separated from the allele number by a dash ("-").
Note: The locus name should exactly match the labels used in the profiles.txt file.

3. Reference file: should be in fasta format and  should be saved as reference.seq. reference.seq should be indexed by using the makeblastdb command:
	makeblastdb -in reference.seq -dbtype nucl  -out reference


Examples are given in the MLST_data directory.

##Usage
-------

	usage: MOST.py [-h]
	version 1-0, date 19/11/2015.

MOST is a modified version of SRST version 1 script (http://sourceforge.net/projects/srst/files/?source=navbar), modification made by  Anthony.Underwood@phe.gov.uk and Rediat.Tewolde@phe.gov.uk.

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

##Running MOST.py

	A. To determine ST
	------------------------------
		Run the script by using the following 2 pathways:
		1.  The script can be run by specifying a workflow and full path to the input directory. The input directory should contain a pair of fastq files.
			e.g MOST.py -w <workflow name>  -i < full path to the input directory>
			Note:If you are using a workflow you must specify an input directory

		2. The script can be run by specifying the full path of the paired fastq files, full path to the directory containing the alleles fasta and profiles.txt files and output directory
			e.g MOST.py -1 <full path to the fastq file> -2 <full path to the fastq file>  -o   <full path to the output directory>  -st <full path to the directory containing the alleles fasta and profiles.txt files>
			Note: If you are passing full fastq paths to the fastq -1 and fastq -2 params, you need also to pass to the -st|--profile_file_directory  param a path to a reference_dir containing
			reference fasta files for each gene of interest and a profiles.txt files.


	B. To determine ST and infer salmonella serotype
	--------------------------------------------------
		Run the script by using the following 2 pathways:
		1.  The script can be run by specifying a workflow, full path to the input directory and set a serotype flag to True. The input directory should contain a pair of fastq files.
			e.g MOST.py -w <workflow name>   -i <path to the input directory> -serotype True
			Note:If using a workflow you must specify an input directory

		2. The script can be run by specifying the full path of the paired fastq files, full path to the directory containing the alleles fasta and profiles.txt files, output directory and set a serotype flag to True.
			e.g MOST.py -1 <full path to the fastq file> -2 <full path to the fastq file> -o  <full path to the output directory>  -st <full path to the directory containing the alleles fasta and profiles.txt files> -serotype True
			Note: If you are passing full fastq paths to the fastq -1 and fastq -2 params, you need also to pass the -st|--profile_file_directory  param, a path to a reference_dir containing
			reference fasta files for each gene of interest and a profiles.txt file.



##Output files
------------

1. Sample_ID_MLST results.csv: describes allele designation and associated confidence quality metrics for each assignment

2. Sample_ID_results.xml: The XML file will report
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

		i.	The "GREEN traffic light" indication is assigned if the
			 max percentage non consensus depth  < 15%  and
			 Complete pileup= TRUE and
			 Minimum consensus depth > 2 and
			 Percentage coverage =100  and
			 ST not "Failed(incomplete locus coverage)

		ii.	The RED traffic light indication is assigned if the
			 Complete pileup= FAIL or
 			 Percentage coverage < 100 or
`			 ST  is  "Failed(incomplete locus coverage)

		iii.	The "AMBER traffic light" indication is assigned if the there is no exact fit which matches either GREEN or RED



##In Galaxy
-----------

MOST is also available as a tool for installation into your own instance of Galaxy (http://galaxyproject.org). Search for 'phemost' in the Galaxy toolshed and install the tool together with the corresponding package and the data manager. Below are detailed instructions on how to do this:

1. Get your own Galaxy server.
    - If you already have a Galaxy server where you can install your own tools from the Galaxy Toolshed continue with point 2.
    - The prerequisites for installing your own Galaxy server are:**
        * A workstation running a Linux operating system
        * Min. 4GB of RAM (better 8GB or more)
        * Root access to that machine.
        * git (optional)
        * Python 2.6 or 2.7
        * GNU Make, gcc to compile and install tool dependencies
    - The last three are standard on most contemporary Linux installations.
    - To get your own Galaxy, please go to the Galaxy Wiki here [https://wiki.galaxyproject.org/Admin/GetGalaxy] and follow the instructions in the sections from "Get the Code" to "Become an Admin" (including).

2. Installing MOST
    - There is a small number of Linux packages that need to be installed on your machine, so that MOST can be installed and run propery. These are development libraries for xml which are required for the installation of the lxml Python library for reading and writing xml files.
    - On some popular Linux distributions the commands to install them are:

    ```
    UBUNTU: sudo apt-get install libxml2-dev libxslt1-dev python-dev
    FEDORA: sudo dnf install libxslt-devel libxml2-devel python-devel
    OPENSUSE: sudo zypper install libxml2-devel libxslt-devel python-devel
    ```

    - Make sure that you can access ftp sites from the command line running your Galaxy server. This is normally enabled by default, but sometimes requires an additional proxy setting. Try this command:

    ```
    wget ftp://ftp.gnu.org/gnu/libtool/libtool-2.4.tar.gz
    ```

    - If that does not download a file named ‘libtool-2.4.tar.gz’ to your current folder, speak to your local systems’ administrator.
    - In your config/galaxy.ini files, set

    ```
    tool_dependency_dir = tool_dependencies
    ```

    - Restart your Galaxy
    - In your Galaxy, click on ‘Admin’ in the main menu at the top.
    - Select ‘Search Tool Shed’ from the menu on the left hand side.
    - Click on the little black triangle (the context menu) next to ‘Galaxy Main Tool Shed’ and select ‘Browse valid repositories’.
    - Type ‘phemost’ [sic] into the search box and press Enter.
    - You should see the “package_phemost_1_0” tool package, the data manager "data_manager_phemost" and the “phemost” tool. Select ‘Preview and install’ from the context menu of the "phemost" tool.
    - Click ‘Install to Galaxy’.
    - Type ‘PHE TOOLS’ into the ‘Add new tool panel section:’ textbox.
    - Click ‘Install’.
    - You will be presented with a long list of packages that need to be installed. This will take a while. Wait until everything is green. If nothing happens for a little while, try reloading the page.

    In your admin tool panel the menu item “Manage installed tools” was added. You can check the status of your installed tools there.

    You also need to install the "data_manager_phemost" for downloading MLST reference data from pubmlst.org for your organism of choice. This is done as above for the "phemost" tool.

3. Use MOST in Galaxy

    First you need to download reference data for your favourite organism, then you can upload your short read data and start tying your samples.

    - Click on 'Admin' in the main menu bar.
    - Click on 'Local data' on the left hand side under 'Data'.
    - Click on 'Fetch MLST Data' under 'Run Data Manager Tools'.
    - Select your organism from the top down menu at the top.
    - Enter an accession number for a complate reference genome for your organsim. Examples for commonly MLST'ed pathogenic bacteria are shown in the table on the same page. The reference genome is automatically downloaded and used to extract the flanking regions of the typing genes.
    - If you organism is not in this table please go to the NCBI Nucleotide archive (http://www.ncbi.nlm.nih.gov/nuccore/) and search for "<species name> complete genome" to find a suitable accession number to use.
    - Click 'Execute'. A 'Fetch MLST Data' job will be added to your history. When the job finished successfully, MOST is ready to use.

    You need to upload your data to a Galaxy history to use it. There are multiple options depending on your local Galaxy configuration. If you have followed the instructions above under ‘Get your own Galaxy’ the only available option is uploading from your local harddrive. When doing this, please make sure your fastq files are in the ‘fastq’ format.

    - Click on the MOST tool under 'PHE TOOLS' on the left hand side.
    - Select your read data in the first and second drop down menu.
    - Your organism should appear in the third drop down menu together with a date stamp of when it was downloaded.
    - Select it and click 'Execute'.
    - A CSV and an XML result file are added to your history.

## Contact
-----------
Anthony.Underwood or Rediat Tewolde
Anthony.Underwood@phe.gov.uk
Rediat.Tewolde@phe.gov.uk
