#PHE MLST software (PMS) is a modified version of SRST version 1 script (http://sourceforge.net/projects/srst/files/?source=navbar), 
#modification made by  Anthony.Underwood@phe.gov.uk and Rediat.Tewolde@phe.gov.uk.
import os
import os.path
import sys
import pickle
import re
from collections import Counter
import log_writer
from utility_functions import *


def get_scores(tmp_dir, profile_file_directory, stderr_log_output):
	
	"""
	Function
	(1) Parse through the pileup file to capture quality metric values

	The option for method:
	tmp_dir[str]: the path to  tmp directory
	profile_file_directory[str]: the path to where  reference.seq, profile.txt and
	the locus variant sequences (*.fas) files located
	
	Return
	results[dict]: key = locus name and  value = quality metric values
	database[dict]: keys correspond to ST numbers and values are array of
	locus variant numbers correspond to ST
	LocusList[list]: list of locus names
	"""
	
	ranges = pickle.load(open(os.path.join(tmp_dir, "ranges.pkl")))
	
	(database, locusList,list_of_all_allele_numbers_tuple) = try_and_except(stderr_log_output,
																			get_profiles,
																			profile_file_directory)
	
	results = try_and_except(stderr_log_output,
							 score,
							 tmp_dir,
							 locusList,
							 ranges)
	
	return results, database, locusList, list_of_all_allele_numbers_tuple


def get_profiles(profile_file_directory):
	
	"""
	Function
	Create a dictionary from profile.txt file in which keys correspond to ST numbers and values
	are array of locus variant numbers correspond to ST
	
	The option of the method
	profile_file_path[str]:  the path to where reference.seq, profile.txt and
	the  locus variant sequences (*.fas) files location
		
	Return
	database[dict]: keys correspond to ST numbers and
	values are array of locus variant numbers correspond to ST
	LocusList[list]: list of locus names
	"""
	
	profile_file_path = profile_file_directory+ "/profiles.txt"
	list_of_all_allele_numbers_tuple = []
	database = None
	locusList = []
	
	for l in open(profile_file_path):
	    if database is None:
		database = {}
		locusList = l.split()[1:]
		continue
	    t = l.split()
	    st = t[0]
	    v = ' '.join([s for s in t[1:]])
	    if v in database:
		print >> sys.stderr, 'sequence type ' + str(st) + ' is a duplicate of ' + str(database[v])
	    database[v] = st
	    covert_string_to_tuple_list_of_allele_numbers = tuple(int(x) for x in re.findall("[0-9]+", v)) 
	    list_of_all_allele_numbers_tuple.append(covert_string_to_tuple_list_of_allele_numbers)
		
	return (database, locusList, list_of_all_allele_numbers_tuple) 


def score(tmp_dir, locusList, ranges):
	
	"""
	Function
	(1) Parse the pileup file
	(2) Calculate quality metric values
	
	
	The option for method:
	tmp_dir[str]: the path to  tmp directory
	LocusList[list]:  list of allele name 
	ranges[str]: start and end position of locus variant sequence
	
	Return
	results[dict]: key = locus base name and  value = quality metric values
	"""
	
	loc = ''
	pos = 1
	count_indel = 0
	holes = 0
	snps = 0
	covMax=combined_covMax=covSum=covSum2= 0 
	covMin = combined_covMin =99999
	percentage_coverages =[]
	snpList = []
	indelList = []
	results = {} 
	
	pileup_file =  os.path.join(tmp_dir, 'all.pileup')
	for l in open(pileup_file):
		t = l.split()
		if loc == '':
			loc = t[0]  
			pos = ranges[loc][0] + 1  
		if t[0] != loc:
			results =GenerateResult(ranges,
									holes, locusList,
									loc,snps,count_indel,
									snpList, indelList,
									percentage_coverages,combined_covMin,
									combined_covMax, covMin, covMax,covSum, results)
			# reset locus vars
			loc = t[0] 
			pos = ranges[loc][0] + 1 
			count_indel = 0
			holes =snps=covMax=combined_covMax=covSum=covSum2= 0 
			covMin =combined_covMin= 99999
	                snpList = []
			indelList = []
			percentage_coverages =[]
		here = int(t[1])
		if here - 1 < ranges[loc][0]:  
			continue
		elif here - 1 >= ranges[loc][1]: 
			continue
		while pos < here: 
			holes += 1 
			pos += 1

		v, indel, array_of_all_indels,most_common_indel = pile(t[2], t[4])
		x = v.items()
		x.sort(lambda a,b: compGreater(t[2], a, b))
		
		if x[0][0] != t[2].lower():
			snps += 1
			snpList.append((pos,t[2],v));
		c = x[0][1] 
		cov= int(most_common_indel)/float(t[3]) 
		if cov > 0.5:  
                    count_indel += 1
                    indel_type = Counter(array_of_all_indels) 
                    indel_type = indel_type.items()
                    indelList.append((int(pos),t[2], indel_type))
		covSum += c 
		covSum2 += c * c
		if c > covMax:
			covMax = c
		if c < covMin:
			covMin = c
		combined_c = x[0][1] + x[1][1] + x[2][1] + x[3][1] 
		if combined_c > combined_covMax:
			combined_covMax = c 
		if combined_c < combined_covMin:
			combined_covMin = c 
		
		n = int(t[3]) 
		js = []
		for (_,j) in x[1:]: 
			js.append(j) 
		percentage_coverage = sum(js)/float(n)*100 
		percentage_coverages.append(round(float(percentage_coverage),2))
		pos = here + 1
	results =GenerateResult(ranges,
							holes,
							locusList,loc,
							snps,count_indel,
							snpList,indelList,
							percentage_coverages,combined_covMin,
							combined_covMax, covMin, covMax,
							covSum, results)
	
	return results


def GenerateResult(ranges, holes, locusList, loc, snps, count_indel, snpList, indelList, percentage_coverages, combined_covMin, combined_covMax, covMin, covMax, covSum, results):
	
	"""
	Function
	If  no gap(holes) in the pileup file report: 
	locus: the locus names(e.g: gki)
	var: variant number (e.g: 2)
	snps:  number of mismatches between the readset and locus
	count_indel: number of mismatches between the readset and locus
	loc: locus variant (e.g: gki-2)
	snpList
	indelList
	quality metric values
	
	"""
	
	nameSep ="-"

	if  len(percentage_coverages) > 0 and holes == 0:
		
		locus_length = int(ranges[loc][1]) - int(ranges[loc][0])
		m = re.search('([^'+nameSep+']+)'+nameSep+'?([0-9]+)', loc) # e.g loc= AROC-4, nameSep : "-"
		locus = m.group(1) 
		if locus not in locusList:
			print "Locus " + locus + " from sequence file not recognised in ST file."
		var = m.group(2) 
		
		max_percentage_of_non_consensus_bases = max(percentage_coverages)
		number_of_time_percentage_coverage_value_calculated =len(percentage_coverages)
		percentage_coverage_number_of_time_percentage_coverage_value_calculated = int(number_of_time_percentage_coverage_value_calculated/float(locus_length)*100)
		
		covStats = (max_percentage_of_non_consensus_bases, combined_covMin, combined_covMax, covMin, covMax, int(100 * covSum / number_of_time_percentage_coverage_value_calculated) / 100.0) #int(100 * covSum / n) / 100.0 is the average
		res = (locus, var, snps, loc, snpList, covStats, indelList, count_indel,number_of_time_percentage_coverage_value_calculated,locus_length,percentage_coverage_number_of_time_percentage_coverage_value_calculated)
		if locus not in results:
			results[locus] = [] 
		results[locus].append(res)
		
	return results


def compGreater(r, a, b):
	

	if a[1] == b[1]: 
		if a[0] == r:
			return -1
		elif b[0] == r:
			return 1
		
	return cmp(b[1], a[1])


def pile(r, s): 
	
	"""
	Function
	Parse through the pileup file and count the number of observed base 
	
	The option for method
	r[str]: Reference base 
	s[str]: depth of coverage
		
	Return
	v[dict]:  number of observed base  at each locus variant position. 
	indel[str]
	array_of_all_indels[list]
	"""
	
	r = r.lower() 
	v = {'a':0, 'c':0, 'g':0, 't':0} 
	s = s.lower()
	most_common_indel = 0
	
	array_of_all_indels = re.findall(r'[\+\-][0-9][0-9]*[A-Za-z]', s)
	if len(array_of_all_indels) > 0:
		most_common_indel = Counter(array_of_all_indels).most_common(1)[0][1]
	indel = len(array_of_all_indels)
	for list_indel in array_of_all_indels:
		number_of_indelsBases = int(list_indel[1:-1])
		position_of_indel = s.find(list_indel)
		s = s[0:position_of_indel] + s[position_of_indel+number_of_indelsBases+len(list_indel)-1:]
	
	skip = False
	for c in s:  
		if skip == True:
			skip = False
		elif c == '.' or c == ',': 
			v[r] += 1 
		elif c == '^': 
			skip = True
		elif c == '$': 
			pass
		elif c == '+' or c == '-':
			pass
		elif c.lower() in v:
			v[c.lower()] += 1
	
	return v, indel, array_of_all_indels, most_common_indel









	

	
