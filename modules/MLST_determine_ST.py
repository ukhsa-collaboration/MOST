#Metric-Oriented Sequence Typer (MOST) software is a modified version of SRST version 1 script (http://sourceforge.net/projects/srst/files/?source=navbar), 
#modification made by  Anthony.Underwood@phe.gov.uk and Rediat.Tewolde@phe.gov.uk.
import sys
import os
import pprint
import log_writer
from utility_functions import *


def DetermineST(stderr_log_output, database, locusList, list_of_all_allele_numbers_tuple, FinalResults):
	
	"""
	Function
	(1) Determine ST values
	(2) Write quality metric values and ST value to MLST_result.csv and result.xml files

	"""
	
	st = try_and_except(stderr_log_output,
						determine_ST,
						database,
						locusList,
						list_of_all_allele_numbers_tuple,
						FinalResults)
	
	return st

	
def determine_ST(database, locusList, list_of_all_allele_numbers_tuple, FinalResults):
	
	"""
	Function
	Determine ST values :
	(1) Determine ST
	(2) If a novel combination of known alleles is detected,the ST is reported as 'NOVEL'
	(3) If one or more loci did not match any existing locus variants, the allele value is falgged by '*',
	the ST is reportd as 'NOVEL_allele', followed by the number of inexact loci
	(4) If one allele variant could not be identified with high confidence, the ST is reported as 'Failed(incomplete locus coverage)'
	
	The option for method:
	database[dict]: keys correspond to ST numbers and values are array of locus variant numbers correspond to ST
	locusList[list]:  List of allele names
	FinalResults
	
	Return
	st[str]: ST value
	"""
	
	st=  None
	locusVariants = [] 
	for locus in locusList:
		locusVariants.append(FinalResults[locus]['VariantNumberHash']) 
	locusVariantsConcat = ' '.join(locusVariants) 
	st = '' 
	closestLocusVariants = []
	if locusVariantsConcat in database:
		st = str(database[locusVariantsConcat]) 
		
	elif '-' in locusVariants: 
		uncertainLoci = locusVariants.count('-') 
		locusVariantsClosest = [] 
		for locus in locusList:
			locusVariantsClosest.append(FinalResults[locus]['ClosestVariantNumber'])
		locusVariantsClosestConcat = ' '.join(locusVariantsClosest)
		
		if '-' in locusVariantsClosest: 
			st = "Failed(incomplete locus coverage)"
		else:
			if locusVariantsClosestConcat in database:
				if int(uncertainLoci) == 1: 
					st =  "NOVEL allele. Closest ST:"+str(database[locusVariantsClosestConcat])  + "(SLV)"
				elif int(uncertainLoci)==2:
					st =  "NOVEL alleles. Closest ST:"+str(database[locusVariantsClosestConcat]) + "(DLV)"
				elif int(uncertainLoci)> 2:
					st =  "NOVEL alleles. Closest ST:"+str(database[locusVariantsClosestConcat]) + "(MLV)"
			else:
				if int(uncertainLoci) == 1:
					st = "NOVEL allele. cannot determine closest ST " + "(SLV)"
				elif int(uncertainLoci)==2:
					st = "NOVEL alleles. cannot determine closest ST " + "(DLV)"
				elif int(uncertainLoci)> 2:
					st = "NOVEL alleles.  cannot determine closest ST" + "(MLV)"
	else:  
		st =  "NOVEL ST"
		if  st == "NOVEL ST":
			allele_variants = locusVariants 
			st = identify_SLV_profile( database, allele_variants,list_of_all_allele_numbers_tuple)	
	
	return st


def identify_SLV_profile(database, allele_variants, list_of_all_allele_numbers_tuple):
	
	
	allele_variants = map(int, allele_variants) 
	list_of_all_common_alleles = []
	locus_with_common_allele_number = []
	closest_match = []
	for list_of_allele in list_of_all_allele_numbers_tuple:
		for variant_number_from_novel_ST, each_variant_number_from_list_of_alleles in zip(allele_variants,list_of_allele):
			if variant_number_from_novel_ST ==  each_variant_number_from_list_of_alleles:
				locus_with_common_allele_number.append(variant_number_from_novel_ST)
		common_alleles = list_of_allele, len(locus_with_common_allele_number), locus_with_common_allele_number
		list_of_all_common_alleles.append(common_alleles)
		locus_with_common_allele_number = []
	single_locus_variables = filter(lambda x: x[1]==6 , list_of_all_common_alleles)
	for single_locus_variable in single_locus_variables:
		clonal_complex= list(single_locus_variable[0])
		allele_numbers = ' '.join(str(MLST_profile) for MLST_profile in clonal_complex)
		ST_value = str(database[allele_numbers])
		ST_value_and_MLST_profile = ST_value
		closest_match.append(ST_value_and_MLST_profile)
	
	if len(closest_match) > 0:
	    st =  "Novel ST. Closest ST:" + closest_match[0] + " (SLV)"
	elif len(closest_match) == 0:
	    st = "NOVEL ST. (no SLV)"
	
	return st 

