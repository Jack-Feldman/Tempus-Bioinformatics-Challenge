#!/bin/env python3

"""
Prototype variant annotation tool created for the Tempus
Bioinformatics Technical Challenge.
Author: Jacob Feldman

Input: VCF file
Output: TSV file containing annotations on each variant. If no output is 
specified program will print table to stdout.

"""

import sys
import requests
import argparse
import os.path

#Add arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="File path to VCF.")
parser.add_argument("-o", type=str, help="Output file name.")
args = parser.parse_args()

#Check if VCF file exists
if not args.i:
	print("No input file specified.")
	exit()

elif not os.path.isfile(args.i):
	print("File specified does not exist.")
	exit()

def get_annotations(hit):
	"""
	Function to retrieve relevant annotations from each variant hit
	Input: Line in VCF file
	Output: Tab separated annotations
	"""
	#Retrieve information about the variant
	variant = hit.split("\t")
	chrom = variant[0]
	pos = variant[1]
	varID = variant[2]
	ref = variant[3]
	alt = variant[4]

	#Create dictionary of info field, value pairs from the INFO section
	info = variant[7].split(";")
	infoDict = dict([x.split('=') for x in info])

	"""Type of variation. Need to list most deleterious if there are more 
	than one type. Here I am unsure the rank order of snp, ins, or del 
	in terms of deleterious possibility. I think it depeds on the specific
	nature of the mutation."""
	typeOfVar = infoDict['TYPE'].split(',')
	if 'complex' in typeOfVar:
		typeOfVar = 'complex'
	else:
		typeOfVar = typeOfVar[0]

	#Depth of sequence coverage at the site of variation
	depth = infoDict["DP"]

	#Number of reads supporting the variant
	numReads = infoDict["AO"]

	#Percentage of reads supporting the variant versus those supporting reference reads
	try: 
		varPercent  = round(float(infoDict["AO"]) / float(infoDict["RO"]) * 100, 4)

	except ZeroDivisionError: #Instance when there are zero reads supporting the reference
		varPercent = "n/a"

	except ValueError: #Instance when there are >1 values for AO or RO. I am unsure why there are more than one count values
		varPercent = "n/a"

	#Allele frequency of variant from Broad Institute ExAC Project API
	#Also inluded is the gene(s) with which the variant is associated with
	res = requests.get(f"http://exac.hms.harvard.edu/rest/variant/variant/{chrom}-{pos}-{ref}-{alt}")
	variantData = res.json()

	try:
		alleleFreq = round(variantData["allele_freq"], 4)
		gene = ','.join(variantData["genes"])

	except KeyError: #Option if there is no data from the api
		alleleFreq = "No data"
		gene = "No data"


	return f"{chrom}\t{pos}\t{varID}\t{typeOfVar}\t{depth}\t{numReads}\t{varPercent}\t{alleleFreq}\t{gene}\n"
#--------------------------------------------------------------------------------------------------------------

#Create basic output table structure
output_table = "Chrom\tPos\tID\tType\tDepth\tNum Reads Supporting Variant\tPercent Variant/Reference Reads\tAllele Freq of Variant\tGenes\n"

#Process VCF file
with open(args.i, "r") as vcf_file:

	line = vcf_file.readline()

	while line.startswith("##"): #Bypass header lines
		line = vcf_file.readline()

	#Write to file if output file is specified, else print to stdout	
	if args.o:
		output = open(args.o, "w")
		output.write(output_table)

		for line in vcf_file: 
			output.write(get_annotations(line))

		output.close()

	else:
		print(output_table)
		for line in vcf_file: 
			print(get_annotations(line))

vcf_file.close()

"""
Thank you Tempus bioinformatics team for taking the time to review my code!
-Jacob
"""