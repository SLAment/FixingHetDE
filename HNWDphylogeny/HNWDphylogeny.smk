## NACHTphylogeny: Phylogeny of the NWD family
#############################################################################
#############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# 2024-03-07 - 2025-04-17
# +++++++++++++++++++++++++++++++++++++++++++++++++

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# -------------------------------------------------
# Configuration file
configfile: "config/config.yaml"
# -------------------------------------------------
# Input files
MASTERalignment = config["MASTERalignment"]
LOCI = config["LOCI"]
# -------------------------------------------------

rule all:
	input:
		expand("results/{locus}.tre", locus = LOCI)

rule break_alignment:
	input:
		fas = MASTERalignment
	output:
		het = "alignments/HET.fa",
		nacht1 = "alignments/NACHT1.fa",
		nacht2 = "alignments/NACHT2.fa",
		nacht3 = "alignments/NACHT3.fa",
	params:
		# Column position for the end of the two first domains (base 1)!!
		## aa2
		NACHTstart = 556, # "RLLWINGDPGKGK" the start
		NACHTend = 1166 # "WISTISVVEAEWN" the end
	run:
		NACHTstart = params.NACHTstart - 1
		NACHTend = params.NACHTend

		# Read the alignment from the FASTA file
		alignment = AlignIO.read(input.fas, "fasta")

		# Rename the sequences in the alignment to keep it simple.
		# Function to rename sequence IDs
		def rename_id(record, delimiter="__"): # All sequences have the format genename__more_info
			new_id = record.id.split(delimiter)[0]
			return SeqRecord(record.seq, id=new_id, description="")

		# Rename all sequences in the alignment
		renamed_records = [rename_id(record) for record in alignment]
		renamed_alignment = MultipleSeqAlignment(renamed_records)

		# Filter out sequences of genes without a HET domain
		genes_sans_HET = ['nwdp-2', "nwd1", "nwd2", "nwd3", 'nwd5', "nwd6"]
		HET_records = [record for record in renamed_alignment if record.id not in genes_sans_HET]
		HET_alignment = MultipleSeqAlignment(HET_records) # Make it a MultipleSeqAlignment object again

		# Calculate the size of each part
		part_size = int(round((NACHTend - NACHTstart) / 3, 0))

		## Create new alignments
		# the first `;,` means "take all rows" (i.e., all sequences in the alignment).
		HETpart = HET_alignment[:, :NACHTstart] # take columns from the start (index 0) up to, but not including, NACHTstart
		part1 = renamed_alignment[:, NACHTstart : NACHTstart+part_size] 
		part2 = renamed_alignment[:, NACHTstart + part_size: NACHTstart + 2*part_size]
		part3 = renamed_alignment[:, NACHTstart + 2*part_size: NACHTend]

		## The HET alingment now has a lot of columns with only gaps
		columns_to_keep = []

		# Iterate over each column
		for i in range(HETpart.get_alignment_length()):
			column = HETpart[:, i]  # Get characters in the i-th column

			# Check if all characters in the column are gaps '-'
			if not all(c == '-' for c in column):
				columns_to_keep.append(i)

		# Create a new alignment with the selected columns
		clean_records = []
		for record in HETpart:
			clean_seq = ''.join(record.seq[i] for i in columns_to_keep)
			clean_records.append(record.__class__(clean_seq, record.id, "", "")) # the last two "" are `name` and `description`
		
		HETpart_clean = MultipleSeqAlignment(clean_records)

		# Save the parts to new FASTA files
		AlignIO.write(HETpart_clean, output.het, "fasta")
		AlignIO.write(part1, output.nacht1, "fasta")
		AlignIO.write(part2, output.nacht2, "fasta")
		AlignIO.write(part3, output.nacht3, "fasta")


rule trimAl:
	input:
		"alignments/{locus}.fa"
	output:
		"trimAl/{locus}_trimAl.fa"
	threads: 1
	shell:
		"trimal -in {input} -out {output} -htmlout trimAl/{wildcards.locus}_trimAl.html -strictplus"

rule IQTreeConcat:
	""" Make a tree of the concatenated alignment """
	input:
		"trimAl/{locus}_trimAl.fa"
	output:
		"IQTree/{locus}.treefile"
	params:
		bootstraps = 100, 
	threads: 8
	shell:
		"iqtree -s {input} -m MFP -seed 1234 -b {params.bootstraps} -nt {threads} -pre 'IQTree/{wildcards.locus}'"
		# -b <#replicates>	 Bootstrap + ML tree + consensus tree (>=100)
		# -bb Ultrafast bootstraps (UFBoots)
		# -bnni to reduce the risk of overestimating branch supports with UFBoot due to severe model violations.

rule get_results:
	input:
		"IQTree/{locus}.treefile"
	output:
		"results/{locus}.tre"
	shell:
		"cp {input} {output}"


