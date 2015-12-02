
import pickle
import sys
import os
import os.path
import pprint
import csv
import glob
from lxml import etree
from collections import Counter
import log_writer
from utility_functions import *


def print_result(stderr_log_output, st, tmp_dir, output_directory, profile_file_directory, ids, workflow_name, version, locusList, FinalResults, infer_salmonella_serotype):

    """
    Function
    Write  quality metric values and ST value to MLST_result.csv and results.xml files

    The option for method:
    st[str]: ST value
    tmp_dir[str]: the path to  tmp directory
    output_directory[str]: the path to the result file
    profile_file_directory[str]: the path to where  reference.seq, profile.txt and Locus variant sequences (*.fas) files located
    workflow_name[str]: tested isolate, e.g. salmonella_typing
    ids[str]: sample unique identifier number
    version[str]: version number
    LocusList[list]:  list of allele names
    """
  
    success, xml_values,FinalResults = try_and_except(stderr_log_output,
                                                      ST_and_CoverageStat,
                                                      st,
                                                      tmp_dir,
                                                      output_directory,
                                                      profile_file_directory,
                                                      workflow_name,
                                                      locusList,
                                                      FinalResults,
                                                      infer_salmonella_serotype)

    if success == False:
        try_and_except(stderr_log_output,
                       report_MLST_result_in_xml_file,
                       xml_values,
                       output_directory,
                       ids,
                       workflow_name,version,infer_salmonella_serotype)


    if success == True:
        try_and_except(stderr_log_output,
                       report_MLST_result_in_csv_file,
                       output_directory,
                       ids,
                       locusList,
                       FinalResults )
        try_and_except(stderr_log_output,
                       report_MLST_result_in_xml_file,
                       xml_values,
                       output_directory,
                       ids,
                       workflow_name,version,infer_salmonella_serotype)

    for root, dirs, files in os.walk(output_directory):
        for currentFile in files:
            exts=('.fasta', '.pkl', '.sam', '.tmp', '.bt2', '.unsortedbam')
            if any(currentFile.lower().endswith(ext) for ext in exts):
                os.remove(os.path.join(root, currentFile))


def ST_and_CoverageStat(st, tmp_dir, output_directory, profile_file_directory, workflow_name, locusList, FinalResults, infer_salmonella_serotype):


  
    name_of_predicted_serotype = predicted_serotype(st,
                                                    profile_file_directory,
                                                    infer_salmonella_serotype)
    FinalResults["predicted_serotype"] =name_of_predicted_serotype

    allele_variants, minimum_consensus_depth, all_minimum_consensus_depth, mean_consensus_depth, all_mean_consensus_depth, max_percentage_non_consensus_base, average_percentage_coverage, all_max_percentage_non_consensus_base = CoverageStat_metrics(FinalResults, locusList)

    if int(average_percentage_coverage) == int(0):
        xml_values = create_xml_file_for_failed_sample(workflow_name)
        return  False, xml_values, FinalResults

    if int(minimum_consensus_depth) <= int(5) or float(max_percentage_non_consensus_base) >= float(15): 
        st = '*'+st
    if int(average_percentage_coverage) < 100:
        st = "Failed(incomplete locus coverage)"

    FinalResults["ST"] =str(st)

    validation_1,validation_2= methodological_validation(tmp_dir,
                                                         locusList,
                                                         FinalResults)

    if st  != "Failed(incomplete locus coverage)"  and validation_1 == 'TRUE'  and int(minimum_consensus_depth) > int(2) and float(max_percentage_non_consensus_base) < float(15) and average_percentage_coverage == 100:
       traffic_light = "GREEN"

    elif validation_1 == 'FALSE' or average_percentage_coverage < 100 or st == "Failed(incomplete locus coverage)" or validation_2 =='FALSE' or int(minimum_consensus_depth) == int(0) :
        traffic_light = "RED"
    else:
        traffic_light = "AMBER"

    xml_values = { 'value': str(st), 'profile': str(allele_variants),
                    'predicted_serotype':str(name_of_predicted_serotype),
                    'QC': {'traffic_light': str(traffic_light),
                    'max_percentage_non_consensus_base': str(max_percentage_non_consensus_base),
                    'max_percentage_non_consensus_base_for_all_loci' : str(all_max_percentage_non_consensus_base),
                    'mean_consensus_depth_for_all_loci':str(all_mean_consensus_depth),
                    'mean_consensus_depth':str(mean_consensus_depth),
                    'minimum_consensus_depth': str(minimum_consensus_depth),
                    'minimum_consensus_depth_for_all_loci':str(all_minimum_consensus_depth) ,
                    'percentage_coverage':str(average_percentage_coverage) ,'complete_pileup':validation_1}}
    
    return  True, xml_values, FinalResults


def CoverageStat_metrics(FinalResults, locusList):


    all_minimum_consensus_depth = []
    all_max_percentage_of_non_consensus_bases = []
    all_mean_consensus_depth = []
    all_allele_variants =[]
    all_loci_percentage_coverage = []
    
    for locus in locusList:
        all_allele_variants.append(FinalResults[locus]["ReportedVariantNumber"])
        all_mean_consensus_depth.append(float(FinalResults[locus]["mean_consensus_depth"])) 
        all_minimum_consensus_depth.append(int(FinalResults[locus]["minimum_consensus_depth"]))
        all_max_percentage_of_non_consensus_bases.append(float(FinalResults[locus]["max_percentage_of_non_consensus_bases"]))
        all_loci_percentage_coverage.append(int(FinalResults[locus]["percentage_coverage"]))

    minimum_consensus_depth = min(all_minimum_consensus_depth)
    all_minimum_consensus_depth = ','.join(str(each_minimum_consensus_depth) for each_minimum_consensus_depth in all_minimum_consensus_depth) 

    mean_consensus_depth = min(all_mean_consensus_depth)
    all_mean_consensus_depth = ', '.join(str(each_mean_consensus_depth) for each_mean_consensus_depth in all_mean_consensus_depth)

    max_percentage_non_consensus_base = max(all_max_percentage_of_non_consensus_bases)
    all_max_percentage_non_consensus_base = ', '.join(str(each_all_max_percentage_of_non_consensus_bases) for each_all_max_percentage_of_non_consensus_bases in all_max_percentage_of_non_consensus_bases)

    average_percentage_coverage = sum(all_loci_percentage_coverage)/int(len(locusList))
    allele_variants = ','.join(all_allele_variants)
    
    return allele_variants, minimum_consensus_depth, all_minimum_consensus_depth, mean_consensus_depth,all_mean_consensus_depth,max_percentage_non_consensus_base,average_percentage_coverage,all_max_percentage_non_consensus_base


def predicted_serotype(st, profile_file_directory, infer_salmonella_serotype):

    """
    Function
    For salmonella samples predicted serotype

    The option for method:
    st[str]: ST value
    profile_file_directory[str]: the path to where  reference.seq, profile.txt and
    the seven Locus variant sequences (*.fas) files located
    workflow_name[str]: tested isolate, e.g. salmonella_typing

    Return
    name_of_predicted_serotype[str]: predicted serotype
    """

    name_of_predicted_serotype = 'none'
    
    if infer_salmonella_serotype:
        predicted_serotype = serotype(st,profile_file_directory)
        print "predicted_serotype", predicted_serotype
      
        name_of_predicted_serotype = predicted_serotype
        if name_of_predicted_serotype == 'no ST-serotype':
            name_of_predicted_serotype = 'no ST-serotype'
        else:
            name_of_predicted_serotype = (sorted(predicted_serotype.items(), key=lambda t: t[1],reverse=True ))
            name_of_predicted_serotype =  str(name_of_predicted_serotype).strip('[]')
    else:
        name_of_predicted_serotype = 'none'
        
    return name_of_predicted_serotype


def serotype(st, profile_file_directory):

    """
    Function
    For salmonella samples, extract the serotype from mlst_lookup.txt file

    The option for method:
    st[str]: ST value
    profile_file_directory[str]: the path to where  reference.seq, profile.txt and
    the seven Locus variant sequences (*.fas) files located

    Return
    predicted_serotype
    """
    
    serotype_lookup_file = open(profile_file_directory + "/mlst_lookup.txt", 'r')
    all_possible_serotypes = []
    
    for line in serotype_lookup_file:
        st_value= line.split('\t')[0]
        serotype= line.split('\t')[1]
        serotype= serotype.strip()
        try:
            if st == st_value:
                all_possible_serotypes.append(serotype)
        except Exception:
            value = 'no ST-serotype'

    predicted_serotype =  Counter(all_possible_serotypes)
    if len(all_possible_serotypes) == 0:
        predicted_serotype = 'no ST-serotype'
        
    return predicted_serotype


def methodological_validation(tmp_dir, locusList, FinalResults):

    """
    Function
    Validation_1: check if seven flanking regions are extracted successfully
    Validation_2: Check if the last sequence in the pileup file (all.pileup)
    is the same as the last sequence in the reference.fa.fai file

    The option for method:
    tmp_dir[str]: the path to  tmp directory

    Return
    validation[str]: returns TRUE if the last sequence in the pileup file (all.pileup) is
    the same as the last sequence in the reference.fa.fai file else return FALSE
    """
    
    reference_file = glob.glob(tmp_dir+'/reference.fa.fai')
    pileup_file = glob.glob(tmp_dir + '/all.pileup')
    last_line_of_the_reference_file ="None"
    last_line_of_the_pileup_file ="None"
    validation_1 = None
    validation_2 = None


    if len(reference_file) >=1:
        fileHandle = open(reference_file[0],"r" )
        lineList = fileHandle.readlines()
        last_line_of_the_reference_file = lineList[-1].split('\t')[0]
        fileHandle.close()
    else:
        last_line_of_the_reference_file ="None"

    if len(pileup_file) >=1:
        fileHandle = open(pileup_file[0],"r" )
        lineList = fileHandle.readlines()
        last_line_of_the_pileup_file = lineList[-1].split('\t')[0]
        fileHandle.close()
    else:
        last_line_of_the_pileup_file ="None"

    if last_line_of_the_pileup_file == last_line_of_the_reference_file:
        validation_1 = "TRUE"
        value = "("+last_line_of_the_pileup_file +":"+last_line_of_the_reference_file+")"
    else:
        validation_1 = "FALSE"
        value = "(" +last_line_of_the_pileup_file + ":" +last_line_of_the_reference_file+")"


    ranges = pickle.load(open(os.path.join(tmp_dir, "ranges.pkl")))
    for locus in locusList:
        numberOfINDELs= FinalResults[locus]["numberOfINDELs"]
    INDELsListsHash= FinalResults[locus]["INDELsListsHash"]
    if int(numberOfINDELs) >= int(1):
        for (pos, ref, TypeOfINDELs) in INDELsListsHash:
            indel_pos_including_flanking_region = pos
            VariantNumberHash = FinalResults[locus]['VariantNumberHash']
            ClosestVariantNumber = FinalResults[locus]['ClosestVariantNumber']
            locus_variant_number = locus+ "-" +ClosestVariantNumber
            indel_pos_including_flanking_region =  ranges[locus_variant_number][1]
            if indel_pos_including_flanking_region == pos:
                validation_2 = "FALSE"
            else:
                validation_2 = "TRUE"

    return validation_1 , validation_2


def report_MLST_result_in_csv_file(output_directory, ids, locusList, FinalResults):


    indel_pos_including_flanking_region = 0
    all_allele_variants =[]
    SNPData = None
    INDELdata = None
    Results_file = output_directory + "/" +ids +"_MLST_result.csv"
    
    with open(Results_file, "wb") as csv_fp:
        dial = csv.excel()
        dial.lineterminator = '\r\n'
        csvWriter= csv.writer(csv_fp, dialect=dial)
        st_value =FinalResults["ST"]
        csvWriter.writerow(["st value:",st_value])
        csvWriter.writerow(["Predicted Serotype",FinalResults["predicted_serotype"]])
        header_row = ["locus name","allele variant","Percentage coverage",
                      "Max percentage of non consensus bases","minimum total depth",
                      "mean consensus depth", "numberOfSNPs","SNPsLists","INDELs","INDELsLists"]

        csvWriter.writerow(header_row)
        for locus in locusList:
            array =[]
            allele_variant =FinalResults[locus]["ReportedVariantNumber"]
            all_allele_variants.append(allele_variant)
            percentage_coverage = FinalResults[locus]["percentage_coverage"]
            max_percentage_of_non_consensus_bases =FinalResults[locus]["max_percentage_of_non_consensus_bases"]
            minimum_total_depth =FinalResults[locus]["minimum_total_depth"]
            mean_consensus_depth= FinalResults[locus]["mean_consensus_depth"]
            SNPsListsHash= FinalResults[locus]["SNPsListsHash"]
            numberOfSNPs = FinalResults[locus]["numberOfSNPs"]
            numberOfINDELs= FinalResults[locus]["numberOfINDELs"]
            INDELsListsHash= FinalResults[locus]["INDELsListsHash"]

            if int(numberOfSNPs) >= int(1):
                for (pos, ref, dist) in SNPsListsHash:
                    ds = 'A(' + str(dist['a']) + ') C(' + str(dist['c']) + ') G(' + str(dist['g']) + ') T(' + str(dist['t']) + ')'
                    pos = pos -100 # remove flanking region
                    SNPData =  "SNP-position:"+str(pos)+ "   reference base:"+ref + "   SNP type:",str(ds)
            else:
                SNPData = None
            if int(numberOfINDELs) >= int(1):
                for (pos, ref, TypeOfINDELs) in INDELsListsHash:
                    ds =','.join(str(dist) for dist in TypeOfINDELs) 
                    pos = pos -100 # remove flanking region
                    INDELdata =  "INDEL-position:"+str(pos)+"  reference base:"+ref +"   INDEL type:"+ str(ds)
            else:
                INDELdata = None

            if int(numberOfINDELs) >= int(1):
                for (pos, ref, TypeOfINDELs) in INDELsListsHash:
                    indel_pos_including_flanking_region = pos

            array.extend((locus, allele_variant, percentage_coverage,
                          max_percentage_of_non_consensus_bases,
                          minimum_total_depth, mean_consensus_depth,
                          numberOfSNPs, SNPData,
                          numberOfINDELs, INDELdata ))
            csvWriter.writerow(array)


def create_xml_file_for_failed_sample(workflow_name):

    if workflow_name == 'salmonella-typing':
        xml_values = { 'value': "Failed(incomplete locus coverage)",'predicted_serotype':'0',  'profile': "0", 'QC' : {'max_percentage_non_consensus_base_for_all_loci' : "0", 'max_percentage_non_consensus_base': "0",'mean_consensus_depth':"0", 'minimum_consensus_depth': "0",  'percentage_coverage': "0",'traffic_light': "RED",'mean_consensus_depth_for_all_loci': "0", 'minimum_consensus_depth_for_all_loci': "0" ,'complete_pileup':"FALSE"}}
    else:
        xml_values = { 'value': "Failed(incomplete locus coverage)",  'profile': "0", 'QC' : {'max_percentage_non_consensus_base_for_all_loci' : "0", 'max_percentage_non_consensus_base': "0",'mean_consensus_depth':"0", 'minimum_consensus_depth': "0", 'percentage_coverage': "0",'traffic_light': "RED",'mean_consensus_depth_for_all_loci':"0", 'minimum_consensus_depth_for_all_loci': "0" ,'complete_pileup':"FALSE"}}
    
    return xml_values


def report_MLST_result_in_xml_file(xml_values, output_directory,ids, workflow_name, version, infer_salmonella_serotype):

    """
    Function
    writes data to results.xml 
    """
    
    xml_log_file = open(output_directory + "/" + ids + ".results.xml", "w")
    
    with open(output_directory + "/" + ids + ".results.xml", "w") as xml_log_file:
        root = etree.Element("ngs_sample", id = ids)
        workflow = etree.SubElement(root, "workflow", value=workflow_name, version = version)
        results = etree.SubElement(root, 'results')
        result = etree.SubElement(results, "result", type="MLST", value = xml_values['value'])

        for types, values in xml_values.items():
            if types == 'profile':
                etree.SubElement(result, "result_data", type=types, value = values)
            if types == "QC":
                for types,value in values.items():
                    etree.SubElement(result, "result_data", type= "QC_" + types, value = value)
            if infer_salmonella_serotype:
                if types == 'predicted_serotype':
                    etree.SubElement(result, "result_data", type=types, value = values)
        print >> xml_log_file, etree.tostring(root, pretty_print=True)
