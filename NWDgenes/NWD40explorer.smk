### HETwd40explorer: Finding and classifying WD40 repeats with high internal conservation in het genes
#############################################################################
#############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# 2023-03-09 - 2025-01-07
# +++++++++++++++++++++++++++++++++++++++++++++++++

# from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio import SeqIO
import re

# -------------------------------------------------
# Configuration file
configfile: "config/config.yaml"
# -------------------------------------------------
# Input files
wdalignment = config["basealignment"]

# Scripts
WDdescriptor = config["descriptor"]
WD40s_PCA = config["WD40s_PCA"]
WD40s_LOGOS = config["WD40s_LOGOS"]

# Variables
AMINOS = config["AMINOS"]
functional_alleles = ["het-D1", "het-D2", "het-E1", "het-E2", "het-E3"]
# -------------------------------------------------


rule all:
	input:
		"results/WD40_LOGOS.png",
		"results/PCA_paper.pdf"

rule WDdescriptor:
	""" Characterize the WD40 repeats """
	# The alignment is special: only the WD40 repeats, no more, no less, and aligned, but without totally empty columns
	input:
		al = wdalignment,
		WDdescriptor = WDdescriptor
	output:
		report = "reports/WDreport_{aminoacids}.txt",
		fasta = "alignments/WDrepeats_{aminoacids}.fa",
		meta = "alignments/WDrepeats_{aminoacids}_metadata.txt",
	shell:
		"aminos=$(echo {wildcards.aminoacids} | sed 's/-/,/g'); "
		"python {input.WDdescriptor} {input.al} -o 'alignments/WDrepeats_{wildcards.aminoacids}' -a $aminos > {output.report}"


# ------

# Process the report to get an R-friendly table of the functional genes
rule WDparser:
	input:
		report = "reports/WDreport_{aminoacids}.txt",
	output:
		report = "reports/WD40_classification_{aminoacids}.txt"
	run:
		ofile = open(output.report, 'w')
		ofile.write("geneid\tallele\tstrain\trepeatn\trepnumber\tposition\n")

		for line in open(input.report, 'r'):
			if 'seqid' in line:
				continue # it's the header
			elif '-----------------------------------' in line: # We reached the end of the report of each allele
				break
			else:
				tab = line.rstrip("\n").split("\t")
				seqid = tab[0]
				wdclass = tab[4].split(',')

				# Extract more data from the name
				seqidlist = seqid.split('_')
				allele = seqidlist[0]
				strain = seqidlist[2].strip('.R10')
				geneid = allele

				# Make a distinction between the gene and the allele
				if 'het-d' in allele or 'het-D' in allele or 'hetD' in allele:
					geneid = 'het-d'
				elif 'het-e' in allele or 'het-E' in allele or 'hetE' in allele:
					geneid = 'het-e'
				elif 'het-r' in allele or 'het-R' in allele or 'hetR' in allele:
					geneid = 'het-r'

				# Fix some problematic seqs
				if 'CU074315' in seqid:
					allele = 'het-e'
					strain = 'PaS'
					geneid = allele
				elif 'L28125' in seqid:
					allele = 'het-E1A'
					geneid = 'het-e'
					if 'fixedlikeEspagne2002' in seqid:
						strain = 'L28125_fix'
					else:
						strain = 'L28125'
				elif 'FJ897789' in seqid: # The last part is missing from the GenBank sequence
					wdclass.append('?')
			
				# Print in long format for R
				if wdclass[0] == '?':
					wdclass = wdclass[1:]

				count = 0
				for repeat in wdclass:
					if repeat != '?':
						repnumber = repeat[1:]
					else:
						repnumber = 0

					if repnumber == "": repnumber = 0

					count += 1
					newline = f"{geneid}\t{allele.replace('het-', '')}\t{strain}\t{repeat}\t{repnumber}\t{count}\n"
					ofile.write(newline)


# Function to translate nucleotide sequence considering gaps and frameshifts
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

rule clean_alignment:
	""" Make per-gene alignments and a report """
	input:
		fasta = "alignments/WDrepeats_{aminoacids}.fa",
		meta = "alignments/WDrepeats_{aminoacids}_metadata.txt",
	output:
		fasta = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean.fa",	
		hete = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_het-e.fa",	
		hetd = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_het-d.fa",	
		hetr = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_het-r.fa",
		hnwd1 = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_hnwd1.fa",	
		hnwd3 = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_hnwd3.fa",	
		nwd1 = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_nwd1.fa",	
		nwd2 = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_nwd2.fa",	
		nwd3 = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_nwd3.fa",	
		nwd5 = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_nwd5.fa",	
		nwd6 = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_nwd6.fa",	
		nwdp2 = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_nwdp-2.fa",
		report = "reports/WDrepeats_{aminoacids}_counts.txt"
	run:
		tabs = [line.rstrip("\n").split("\t") for line in open(input.meta, 'r')] 			# Read tab file into a list

		# Filter for the sequences that are part of the paper
		accepted_strains = ["Podan2", "PaYp", "PaZp.R10", "PaTgp", "PaWa21m", "PaWa28m", "PaWa46p", "PaWa53m", "PaWa63p", "PaWa87p", "PaWa100p", "PaWa137m"]
		backcrossed_strains = ["CmEmm.R10", "CoEcp.R10", "CoEfp.R10", "ChEhDap.R10", "CaDam.R10", "CsDfp.R10"]

		ans_tab = []
		genedic = {}
		for tab in tabs[1:]: # ignore the header (first line)
			SequenceName = tab[0]
			Strain = SequenceName.split('_')[3]
			Species = tab[3]
			Gene = tab[4]
			
			if Strain in accepted_strains + backcrossed_strains:
				if Strain in backcrossed_strains:
					Allele = tab[5]
					if Allele in functional_alleles:
						status = True # keep track of how acceptable this sequence is
					else:
						status = False
				else:
					status = True

				if status: # If it's ok, save it
					ans_tab.append(SequenceName)

					if Gene in genedic.keys():
						genedic[Gene].append(SequenceName)
					else:
						genedic[Gene] = [SequenceName]

				status = False # reset 

		# This stores in memory
		records_dict = SeqIO.to_dict(SeqIO.parse(input.fasta, "fasta"))

		# All sequences
		total_clean = 0
		total_clean_filtered = 0
		outputfasta = open(output.fasta, "w")
		for seq in ans_tab:
			total_clean += 1
			protein_seq = custom_translate(records_dict[seq].seq)

			# Accept all repeats
			SeqIO.write(records_dict[seq], outputfasta, "fasta")
			total_clean_filtered +=1

			# if 'X' not in protein_seq: # Keep only repeats that are intact
			# 	SeqIO.write(records_dict[seq], outputfasta, "fasta")
			# 	total_clean_filtered +=1

		# print a little report
		reportfile = open(output.report, 'w')
		reportfile.write(f"Total number of processed repeats: {len(records_dict.keys())}\n")
		reportfile.write(f"Total number of clean repeats: {total_clean}\n")
		reportfile.write(f"Total number of clean filtered repeats: {total_clean_filtered}\n")

		# The gene specific alignments
		for Gene in genedic.keys():
			total_gene = 0
			total_gene_filtered = 0

			genefasta = output.fasta.replace('.fa','') + '_' + Gene + '.fa'
			outputfasta = open(genefasta, "w")
			for seq in genedic[Gene]:
				total_gene += 1
				protein_seq = custom_translate(records_dict[seq].seq)
				
				# Accept all repeats
				SeqIO.write(records_dict[seq], outputfasta, "fasta")
				total_gene_filtered += 1

				# if 'X' not in protein_seq: # Keep only repeats that are intact
				# 	SeqIO.write(records_dict[seq], outputfasta, "fasta")
				# 	total_gene_filtered += 1

			reportfile.write(f"Repeat number of {Gene}: {total_gene}\n")
			reportfile.write(f"Repeat number of {Gene} after filtering: {total_gene_filtered}\n")

rule plotPCA:
	""" Make PCA plots of the WD repeats based on paralog identity """
	input:
		meta = expand("alignments/WDrepeats_{aminoacids}_metadata.txt", aminoacids = AMINOS),
		fasta = expand("alignments/PerSpecies/WDrepeats_{aminoacids}_clean.fa", aminoacids = AMINOS)
	output:
		pcaall = "results/PCA_allfam.pdf",
		pcahete = "results/PCA_heteclade.pdf",
		pcahetd = "results/PCA_hetdclade.png",
		paper = "results/PCA_paper.pdf"
	conda:
		"envs/logopca.yaml"
	script:
		WD40s_PCA

rule nc2aa:
	""" I want to have the repeats in amino acids as well """
	input:
		fasta = "alignments/PerSpecies/WDrepeats_{aminoacids}_clean_{gene}.fa"
	output:
		fasta = "alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_{gene}_aa.fa"
	run:
		# Initialize a list to hold the protein sequences
		protein_records = []

		ofile = open(output.fasta, 'w')
		# Read and translate the nucleotide sequences
		for record in SeqIO.parse(input.fasta, "fasta"):
			# Translate the nucleotide sequence to a protein sequence with custom function to tolerate frame-shifts
			protein_seq = custom_translate(record.seq)
			protseq = ">" + record.id + "\n" + protein_seq + "\n" # Using f-strings makes weird errors with appending extra spaces 
			ofile.write(protseq)

rule makeLOGOS:
	input:
		hete = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_het-e_aa.fa", aminoacids = AMINOS),
		hetd = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_het-d_aa.fa", aminoacids = AMINOS),
		hetr = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_het-r_aa.fa", aminoacids = AMINOS),
		hnwd1 = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_hnwd1_aa.fa", aminoacids = AMINOS),
		hnwd3 = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_hnwd3_aa.fa", aminoacids = AMINOS),
		nwd1 = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_nwd1_aa.fa", aminoacids = AMINOS),
		nwd2 = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_nwd2_aa.fa", aminoacids = AMINOS),
		nwd3 = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_nwd3_aa.fa", aminoacids = AMINOS),
		nwd5 = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_nwd5_aa.fa", aminoacids = AMINOS),
		nwd6 = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_nwd6_aa.fa", aminoacids = AMINOS),
		nwdp2 = expand("alignments/PerSpecies/Protein/WDrepeats_{aminoacids}_clean_nwdp-2_aa.fa", aminoacids = AMINOS),
	output:
		eachgene = "results/WD40_LOGOS.png"
	conda:
		"envs/logopca.yaml"
	script:
		WD40s_LOGOS



