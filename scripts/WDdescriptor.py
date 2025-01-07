#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function

# ================== WDdescriptor =================
# Script to characterize the WD domain of the HNWD genes in the Podospora species complex. 
# Inspired on HetRdescriptor.py

# The input is a special alignment of nucleotides but aligned to preserve the
# coding frame of the WD repeats.

# ==================================================
# Sandra Lorena Ament-Velasquez
# 2023/02/22
# +++++++++++++++++++++++++++++++++++++++++++++++++
# 2024/03/19 - v1.32, fixed an issue with the repeat naming for unknown repeats ("?")
# 2024/03/19 - v1.31, finally removed code from original Chevanne et al. 2010 classification
# 2024/02/28 - v1.3, added recognition of nwd1
# ------------------------------------------------------
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

import sys # For reading the input
import re
import os # For the input name
from collections import defaultdict  # This tells python that the dictionary contains a list so you can freely append things to the value of each key
import argparse # For the fancy options
# ------------------------------------------------------

version = 2.0
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Characterize the WD40 diversity of HNWD paralogs *", epilog="Warning: the script doesn't work properly if the input sequences have tracks of IUPAC symbols at the beginning and end of the sequences, or '?' symbols anywhere.") # Create the object using class argparse

# Mandatory options
parser.add_argument('fasta', help="Fasta sequence of just the WD40 motif (the repeats) of an HWND paralog")
# Extra options
parser.add_argument("--output", "-o", help="Base of output name (nothing by default, so reporting to working directory)", default="")
parser.add_argument("--aminoacids", "-a", help="String of comma-separated amino acid residues to classify the WD40 repeats. Deault: 5,6,7,9,25,27", default = "5,6,7,9,25,27", type=str)

# Input from console
try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	fastaopen = open(args.fasta, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()
# except:
# 	print("Usage: python " + sys.argv[0] + " hnwd_master_onlyWDdomain_noemptycols_noGuides_hetsfirst.fa")
# 	print("Version", versiondisplay)
# 	sys.exit(1)

# ------------------------------------------------------

# ---------------------------------
# Get sequences in file
# ---------------------------------
# Function to translate nucleotide sequence considering gaps and frameshifts
# Thanks ChatGPT
# See also https://chatgpt.com/c/eba0aff2-cb16-4047-b1e1-08cbef511c39
def custom_translate(nuc_seq): 
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    protein_seq = []

    # Translate each codon
    for i in range(0, len(nuc_seq), 3):
        codon = str(nuc_seq[i:i+3])
        if len(codon) < 3 or '-' in codon:  # Incomplete codon or indel
            protein_seq.append('X')  # 'X' for unknown amino acid
        else:
            aa = standard_table.forward_table.get(codon, 'X')  # Translate or 'X' if invalid codon
            protein_seq.append(aa)

    return ''.join(protein_seq)

seqdic = {}

# Older version
for seq_record in SeqIO.parse(fastaopen, "fasta"):
	geneseq = seq_record

	# Remove the spaces I have for me to see the WDs better, and the actual
	# frame shifts make them ambiguous bases so that the translation works 
	# Indels are also removed
	ncseq = str(geneseq.seq).upper().strip("-").replace("---","")

	# Make a new record
	geneseq.seq = Seq(ncseq)

	# Finally translate
	aaseq = custom_translate(geneseq.seq)

	# Add to the working dictionary
	seqdic[seq_record.id] = (ncseq, aaseq)

# for seq in seqdic.keys(): print(seq, seqdic[seq])
# print(seqdic)

# ---------------------------------
# Functions
# ---------------------------------
seqidregex = re.compile(r"([h]?nwd[a-zA-Z0-9\-]{0,5}|het-.{1,2})(_{1,2})([a-zA-Z0-9\.]+)_(.*)")

def simplifynames(name):
	namematch = seqidregex.search(name.replace("New|", ""))
	if namematch:
		newname = f"{namematch.group(1).rstrip('_')}__{namematch.group(3)}"
	else:
		newname = name
	return(newname)

# Variables to keep track of numbers and names
WDid_dic = {}

WDcount4names_R = 0
WDcount4names_D = 0
WDcount4names_E = 0
WDcount4names_hnwd1 = 0
WDcount4names_hnwd3 = 0
WDcount4names_nwd1 = 0
WDcount4names_nwd2 = 0
WDcount4names_nwd3 = 0
WDcount4names_nwd4 = 0
WDcount4names_nwd5 = 0
WDcount4names_nwd6 = 0
WDcount4names_nwdp2 = 0
WDcount4names = 0


# # function of user J.F. Sebastian in https://stackoverflow.com/questions/592746/how-can-you-print-a-variable-name-in-python
# # This fails if there is more than one variable name referring to the same value
# def namestr(obj):
#     return [name for name in globals() if globals()[name] is obj]

def assignWDid(WDrepeat, baseseqid):
	global WDid_dic
	global WDcount4names_R
	global WDcount4names_D
	global WDcount4names_E
	global WDcount4names_hnwd1
	global WDcount4names_hnwd3
	global WDcount4names_nwd1
	global WDcount4names_nwd2
	global WDcount4names_nwd3
	global WDcount4names_nwd4
	global WDcount4names_nwd5
	global WDcount4names_nwd6
	global WDcount4names_nwdp2
	global WDcount4names

	aminoacids = args.aminoacids.split(",")
	
	# Substract the chosen amino acids from the WD40 repeat in a specific string
	specificstring = ''
	for amino in args.aminoacids.split(","):
		# print(int(amino), file=sys.stderr)
		specificstring += WDrepeat[int(amino) - 1]

	# WDid = "?"

	if specificstring not in WDid_dic.keys():
		if "het-r" in baseseqid or 'het-R' in baseseqid or 'hetR' in baseseqid:
			WDcount4names_R += 1
			WDid_dic[specificstring] = "r" + str(WDcount4names_R)
		elif "het-e" in baseseqid or 'het-E' in baseseqid or 'hetE' in baseseqid:
			WDcount4names_E += 1
			WDid_dic[specificstring] = "e" + str(WDcount4names_E)
		elif "het-d" in baseseqid or 'het-D' in baseseqid or 'hetD' in baseseqid:
			WDcount4names_D += 1
			WDid_dic[specificstring] = "d" + str(WDcount4names_D)
		elif "hnwd1" in baseseqid:
			WDcount4names_hnwd1 += 1
			WDid_dic[specificstring] = "w" + str(WDcount4names_hnwd1)
		elif "hnwd3" in baseseqid:
			WDcount4names_hnwd3 += 1
			WDid_dic[specificstring] = "y" + str(WDcount4names_hnwd3)
		elif "nwd1" in baseseqid:
			WDcount4names_nwd1 += 1
			WDid_dic[specificstring] = "n" + str(WDcount4names_nwd1)
		elif "nwd2" in baseseqid:
			WDcount4names_nwd2 += 1
			WDid_dic[specificstring] = "m" + str(WDcount4names_nwd2)
		elif "nwd3" in baseseqid:
			WDcount4names_nwd3 += 1
			WDid_dic[specificstring] = "o" + str(WDcount4names_nwd3)
		elif "nwd4" in baseseqid:
			WDcount4names_nwd4 += 1
			WDid_dic[specificstring] = "p" + str(WDcount4names_nwd4)
		elif "nwd5" in baseseqid:
			WDcount4names_nwd5 += 1
			WDid_dic[specificstring] = "q" + str(WDcount4names_nwd5)
		elif "nwd6" in baseseqid:
			WDcount4names_nwd6 += 1
			WDid_dic[specificstring] = "t" + str(WDcount4names_nwd6)
		elif "nwdp2" in baseseqid or "nwdp-2" in baseseqid:
			WDcount4names_nwdp2 += 1
			WDid_dic[specificstring] = "u" + str(WDcount4names_nwdp2)
		else:
			WDcount4names += 1
			WDid_dic[specificstring] = "x" + str(WDcount4names)

	return(WDid_dic[specificstring])

def assigngene(baseseqid):
	allele = baseseqid
	
	generegex = re.compile("(het-|nwd|hnwd|nwdp-)([a-zA-Z0-9]{1,2})([a-zA-Z]?)__") # re.compile("([a-zA-Z0-9-]+)__")
	namematch = generegex.search(baseseqid)

	if namematch:
		allele = namematch.group(1) + namematch.group(2)
	else:
		allele = baseseqid
		# print(allele, file=sys.stderr)
		

	baseseqid = baseseqid.lower() # Make them lowercase to remove confusion
	# The order matters, otherwise it overlaps with nwd1
	genelist = ["het-r", "het-e", "het-d", "hnwd1", "hnwd3", "nwd1", "nwd2", "nwd3", "nwd4", "nwd5", "nwd6", "nwdp2", "nwdp-2"]

	for thisgene in genelist:
		if thisgene in baseseqid:
			gene = thisgene
			if gene == "nwdp2":
				gene = "nwdp-2"
			break
		elif namematch:
			gene = namematch.group(1)
		else:
			gene = "?"
			
	return([gene, allele])

def assignspp(strainID):
	species = "?"
	if "Pa" in strainID:
		species = "anserina"
	elif "Pc" in strainID:
		species = "comata"
	elif strainID in ["CmEmm.R10", "CoEcp.R10", "CoEcm.R10", "CoEfp.R10", "ChEhDap.R10", "ChEhDa.R10", "CaDam.R10", "CsDfp.R10", "CsDfm.R10"]:
		species = "anserina"
	elif "CBS333.63" in strainID or "CBS451.62" in strainID or "CBS237.71" in strainID:
		species = "pauciseta"
	elif "CBS411.78" in strainID:
		species = "pseudoanserina"
	elif "CBS415.72" in strainID:
		species = "pseudocomata"
	elif "CBS112042" in strainID:
		species = "bellae-mahoneyi"
	elif "CBS124.78" in strainID or "CBS253.71" in strainID:
		species = "pseudopauciseta"
	return(species)


# ---------------------------------
# Find the WDs
# ---------------------------------

# Output files
outputfasta = open( f"{args.output}.fa", "w")
outputmetadata = open(f"{args.output}_metadata.txt", "w")

outputmetadata.write("Sequence\tRepeatType\tStrain\tSpecies\tGene\tAllele\tSeqID\tPosition\tEnd\n") # Print a header for the metadata

# REGEX to find WD40 repeats
# StrictWDmotif = re.compile(r"(L|F|I)(E|K|A|Q)([\w]{36})(C|Y|G|R|K|E|H|S)(T|I|F|L|M|V|Q|R|W)(Q|R|H|K|L|E)(T|A|M|K)") # Original, based on Paoletti et al. (2007) 10.1371/journal.pone.0000283
StrictWDmotif = re.compile(r"([\w]{5})(L|F|I)(E|K|A|Q)(G|S)(H|Y)([\w]{31})(S|T)(G)") # WD repeat following Hu et al. (2017) 10.1038/s41598-017-11115-1
# StrictWDmotif = re.compile(r"(G|E|D|L|S|C)(H|Y)([\w]{21})(D|G)([\w]{5})(W|P|R)([\w]{4})(G|E|S|C)([\w]{7})") # WD40 blade following Hu et al. (2017) 10.1038/s41598-017-11115-1

strainregex = re.compile(r"(Pa|Pc|CBS)([\w\.]+)")
strainR10regex = re.compile(r"([a-zA-Z]*)(p|m)(.R10)")
knowngeneregex = re.compile(r'het-[a-zA-Z][0-9]')

# Keep track of the origin of everything
generalcount = 0

# Header of results
sys.stdout.write(f'seqid\tWD40domain_len_aa\tExpected_No_repeats\tCount_Good_WD40\tWD40_classification\tCterm\n') 

for seqid in seqdic.keys():
	baseseqid = simplifynames(seqid)
	REPEATseq = seqdic[seqid][1]

	# Define the boundaries of the WD repeat region
	ExpRep = float(len(REPEATseq))/42 # If there are perfect repeats, this should be a whole number

	# Loop through the expected length of the WD40 C-term:
	chippingREPEATseq = REPEATseq
	chippingREPEATseqNC = seqdic[seqid][0]

	countREPEAT = 0
	REPEATstring = ''
	badREPEATstring = ''
	RPdomainStart = 0 
	for n in range(0, len(REPEATseq)): # Loop the whole sequence (a bit slow if there are sections that are not repeats)
	# for n in range(0, int(ExpRep)): # Loop through all the WD repeats
		WDmatch = StrictWDmotif.search(chippingREPEATseq) # Look for a specific shape of WD repeat as defined above
		if WDmatch:
			currentWD = WDmatch.group(0)
			wdstart, wdend = WDmatch.span()
			
			if 'X' not in currentWD:
				if wdstart > 0: # something is wrong with that repeat
					if len(currentWD) == 42: # It might still be okay but preceed by a broken repeat
						wdid = assignWDid(currentWD, seqid)
						if REPEATstring == '': # This is the first encountered repeat
							REPEATstring += wdid + ','
							RPdomainStart = wdstart # There is seq before the repeat domain in this alignment
						else:				
							REPEATstring += "?," + wdid + ','
					else: # Probably the repeat is broken
						REPEATstring += "?,"
						wdid = "?"
				else: # The Repeat is well-behaved
					wdid = assignWDid(currentWD, seqid)
					REPEATstring += wdid + ','
			else: # There is a stop codon so we can't tedermine what repeat this is
				REPEATstring += "!,"
				wdid = "!"

			# Retrieve sequence into the fasta file
			generalcount += 1 # To make every single repeat unique
			currentWDnc = chippingREPEATseqNC[wdstart * 3 : wdend * 3] # KEY assumption: the base nucleotide alignment is made in such a way that it conserves the triplets
			outputfasta.write(f">{wdid}_{baseseqid}_{generalcount}\n")
			outputfasta.write(f"{currentWDnc}\n")

			# Get the metadata
			strainID = "" # Re-start so it doesn't carry from previous loop
			species = ""

			## Infer the strain and the species
			if "Podan" in baseseqid or "CAL30215.1" in baseseqid:
				strainID = "PaSp"
				species = "anserina"
			elif "PODCO" in baseseqid:
				strainID = "PcTdp"
				species = "comata"
			elif "FJ897789" in baseseqid:
				strainID = "PaA"
				species = "anserina"
			else:
				strainmatch = strainregex.search(baseseqid)
				strainR10match = strainR10regex.search(baseseqid)
				if strainR10match:
					strainID = strainR10match.group(1) + strainR10match.group(2)
					species = "anserina"
				elif strainmatch:
					strainID = strainmatch.group(1) + strainmatch.group(2)
					species = assignspp(strainID)
				else:
					species = "Podospora"

				if knowngeneregex.search(baseseqid): #Some genes have known alleles and those are only within anserina
					species = "anserina"
					if strainID == "": strainID = "Pa"

			outputmetadata.write(f"{wdid}_{baseseqid}_{generalcount}\t{wdid}\t{strainID}\t{species}\t{assigngene(baseseqid)[0]}\t{assigngene(baseseqid)[1]}\t{generalcount}\t{n+1}\n")

			# Move to the next section of the sequence
			chippingREPEATseq = chippingREPEATseq[wdend:] # Move to the next one
			chippingREPEATseqNC = chippingREPEATseqNC[wdend * 3:] # Move to the next one

			countREPEAT += 1 # Count this repeat as whole, for reporting
			currentWD = "" # Re-start so it doesn't carry from previous loop
			strainID = "" # Re-start so it doesn't carry from previous loop
		# else:
		# 	print(chippingREPEATseq)


	# Put something if the REPEATstring is empty
	if REPEATstring == '': REPEATstring = '.'

	# Report
	lenREPEATdomain =  len(REPEATseq) - RPdomainStart - len(chippingREPEATseq)
	resultline = seqid + '\t' + str(lenREPEATdomain) + '\t' + format(lenREPEATdomain/42, '0.2f') + '\t' + str(countREPEAT) + '\t' + REPEATstring.rstrip(',') + '\t' + chippingREPEATseq + '\n'

	sys.stdout.write(resultline) 

outputfasta.close()
outputmetadata.close()

# Reporting
print("-----------------------------------")
print(f"Amino acids used for classification: {args.aminoacids}")
print("Repeat number of het-d:", WDcount4names_D)
print("Repeat number of het-e:", WDcount4names_E)
print("Repeat number of het-r:", WDcount4names_R)
print("Repeat number of hnwd1:", WDcount4names_hnwd1)
print("Repeat number of hnwd3:", WDcount4names_hnwd3)
print("Repeat number of nwd1:", WDcount4names_nwd1)
print("Repeat number of nwd2:", WDcount4names_nwd2)
print("Repeat number of nwd3:", WDcount4names_nwd3)
print("Repeat number of nwd4:", WDcount4names_nwd4)
print("Repeat number of nwd5:", WDcount4names_nwd5)
print("Repeat number of nwd6:", WDcount4names_nwd6)
print("Repeat number of nwdp2:", WDcount4names_nwdp2)
print("Repeat number of Unknown:", WDcount4names)
print("Total", len(list(WDid_dic.values())), '\n')

print(WDid_dic)
