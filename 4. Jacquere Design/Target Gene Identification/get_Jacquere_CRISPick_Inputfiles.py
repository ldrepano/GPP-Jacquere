'''
Laura Drepanos
January 2025
Genetic Perturbation Platform R&D, Broad Institute of MIT and Harvard

This script generates files required to design Jacquere, a CRISPRko Cas9 library targeting all protein coding genes identified by Ensembl, RefSeq, or CHESS 

------------------------------------------------

Jacquere is designed in four CRISPick runs:
1. Targets set of all protein-coding genes with Ensembl as reference genome to target all recognized by GENCODE (v47)
2. Targets set of all protein-coding genes that are recognized by RefSeq (v2024_08) but not GENCODE, using NCBI as reference genome (so that gene symbols are recognized)
3. Targets genes recognized uniquely by CHESS (v3.1.3), using Ensembl as reference genome and supplying the coordinates of the primary CHESS transcripts 
4. Generate 900 intergenic and 100 non-targeting controls with Ensembl as the reference genome

Steps to collect input for CRISPick runs (updated 1/29/2025)
1. Get RefSeq, GENCODE genes (repeat with NCBI, Ensembl the CRISPick reference genome)
	a. Run CRISPick and select all protein coding genes. Since CRISPick is kept up to date, this will yield the most current annotations
	b. Open the sgrna-designs output file in excel and select the “Target Gene Symbol” column

	Steps c.- e. are only required for the Ensembl run (NCBI run does not yield blanks)
	c. You need to remove the blanks in this column for recognition by the HUGO web interface. Under the “Data” excel tab, select “Remove Duplicates”
	d. Command + shift + down-arrow to the first instance of an empty cell, and copy all cells below that up to remove that blank 
	e. Copy the entire column now that all blanks have been removed 

2. Get CHESS genes
	a. Download GRCh38 annotations that map to the primary assembly (fourth Content) from CHESS downloads page (https://ccb.jhu.edu/chess/)
	b. Extract gene names with the following code in command line:

grep -o 'gene_id "CHS.[0-9]*"; gene_name "[a-zA-Z0-9-]*"; gene_type "protein_coding"' chess3.1.3.GRCh38.assembly.gtf | uniq > chess3.1.3_genes_temp.txt
grep -o "[A-Z][a-zA-Z0-9-]*" chess3.1.3_genes_temp.txt | sort >chess3.1.3_genes.txt
	
	Note to not delete the chess3.1.3_genes_temp.txt file: it is used by this script for Step 4e

3. Paste gene name lists into HUGO (https://www.google.com/url?q=https://www.genenames.org/tools/multi-symbol-checker/&sa=D&source=docs&ust=1738188592203904&usg=AOvVaw3ntWOEY22adlG0o3q0zGSB) to get the approved symbol when possible, output as csv

4. Use this python script to 
	a. Take HUGO outputs, replace empty values in approved symbol column with input 
	b. Find overlap of updated approved symbol columns
	c. Get RefSeq gene names not in Ensembl
	d. Get genes unique to CHESS (not present in RefSeq or Ensembl)
	e. Get coordinates of (primary transcript of) unique CHESS genes into input format for CRISPick 

------------------------------------------------

Additional input files for this script: 
"chess3.1.3.GRCh38.assembly.bed": Retrieves the coordinates of CHESS unique genes 
Downloaded CHESS GRCh38 annotations in big bed annotation from github. Converted to bed format with UCSC bigBedtoBed tool.


Outputs:
"refseq_2024_08_not_ensembl_133_genes_proteincoding_assembly.csv": Input for CRISPick run #2
"CHESS3.1.3_unique_protein_coding_genes_coordinates.csv": Input for CRISPick run #3 
"CHESS3.1.3_unique_protein_coding_genes_coordinates_with_identifier.csv": Associates CHESS gene coordinates with identifier. Useful for downstream analysis

'''

import pandas as pd 
import numpy as np
from matplotlib_venn import venn3  
from matplotlib_venn.layout.venn3 import DefaultLayoutAlgorithm
from matplotlib import pyplot as plt 
import gpplot as gpp

gpp.set_aesthetics()


def get_unique_chess_genes():
	#Outputs from Step 3 detailed in file description
	chess_hugo_filepath="Target-Gene-Identification-Data/hgnc-symbol-check_chess3.1.3_GRCh38_assembly_proteincoding.csv"
	ensembl_hugo_filepath="Target-Gene-Identification-Data/hgnc-symbol-check-ensembl113_GRCh38_assembly_proteincoding.csv"
	ncbi_hugo_filepath="Target-Gene-Identification-Data/hgnc-symbol-check_NCBI_2024_assembly_proteincoding.csv"

	chess_hugo=pd.read_csv(chess_hugo_filepath,skiprows=[0])
	ensembl_hugo=pd.read_csv(ensembl_hugo_filepath,skiprows=[0])
	ncbi_hugo=pd.read_csv(ncbi_hugo_filepath,skiprows=[0])

	#Refer to gene by Approved Hugo symbol, unless one does not exist 
	chess_hugo["Resolved gene symbol"]= np.where(chess_hugo["Match type"]!="Unmatched", chess_hugo["Approved symbol"],chess_hugo["Input"])
	chess_genes=chess_hugo["Resolved gene symbol"].tolist()
	ensembl_hugo["Resolved gene symbol"]= np.where(ensembl_hugo["Match type"]!="Unmatched", ensembl_hugo["Approved symbol"],ensembl_hugo["Input"])
	ensembl_genes=ensembl_hugo["Resolved gene symbol"].tolist()
	ncbi_hugo["Resolved gene symbol"]= np.where(ncbi_hugo["Match type"]!="Unmatched", ncbi_hugo["Approved symbol"],ncbi_hugo["Input"])
	ncbi_genes=ncbi_hugo["Resolved gene symbol"].tolist()

	#Generate venn diagram to visualize overlap of genes identified by different catalogs
	v=venn3([set(chess_genes),set(ensembl_genes),set(ncbi_genes)],set_labels=('CHESS 3.1.3','GENCODE47','RefSeq (08-2024)'),
			layout_algorithm=DefaultLayoutAlgorithm(fixed_subset_sizes=(1,1,1,1,1,1,1)),
          	set_colors= ("y","g","b")) 
	#gpp.savefig("11_2024_CHESS_GENCODE_RefSeq_proteincoding_assembly_geneoverlap.pdf")

	#Get list of refseq genes not in ensembl (by their refseq, not hugo name)
	ncbi_unique_genes= [gene for gene in ncbi_genes if gene not in ensembl_genes]
	ncbi_unique_genes=ncbi_hugo[ncbi_hugo["Resolved gene symbol"].isin(ncbi_unique_genes)]["Input"]
	#Generate input file for CRISPick run #2 
	ncbi_unique_genes.to_csv("Target-Gene-Identification-Data/refseq_2024_08_not_ensembl_133_genes_proteincoding_assembly.csv",index=False,header=False)

	#Get list of genes unique to CHESS (by their chess, not hugo name)
	chess_unique_genes= [gene for gene in chess_genes if gene not in ncbi_genes and gene not in ensembl_genes]
	chess_unique_genes=chess_hugo[chess_hugo["Resolved gene symbol"].isin(chess_unique_genes)]["Input"].tolist()
	return chess_unique_genes


def get_chess_coordinates(CHESS3_unique_genes):
	#Get coordinates for transcripts targeting genes uniquely present in CHESS

	'''
	
	
	Notes on converting bed to desired output: 

	Column 12 (0 index) of each row (a transcript) in the bed file indicates the gene that the transcript targets
	Column 3 indicates the entire transcript identifier (CHS.geneID.transcriptID)
	Column 9 indicates the number of exons for that transcript 
	Column 10  indicates block sizes of all exons present in the transcript. 
	Column 11 provides the start positions (relative to the start position in column 1) 
	Column 0 provides the chromosome. 
	CRISPick suggests input format for manual coordinates as NC_000001.11:-:12000-13000;15000-16000 for example, 
	but can resolve the chromosome identifier if you merely input “chr1” for example.   

	'''
	#get CHESS identifiers of CHESS genes
	chess_gene_to_identifier= pd.read_table("Target-Gene-Identification-Data/chess3.1.3_genes_temp.txt",delimiter=" ",header=None)
	chess_gene_to_identifier["identifier"]=chess_gene_to_identifier[1].str[:-1]
	chess_gene_to_identifier["gene_name"]=chess_gene_to_identifier[3].str[:-1]
	CHESS3_unique_gene_identifiers=chess_gene_to_identifier[chess_gene_to_identifier["gene_name"].isin(CHESS3_unique_genes)]["identifier"].tolist()
	
	#read in CHESS3.1.3 bedfile mapped to GRCh38
	chess_bed_file= pd.read_csv("Target-Gene-Identification-Data/chess3.1.3.GRCh38.assembly.bed",delimiter="\t",header=None)
	#subset to only include transcripts targeting protein coding genes unique to CHESS
	unique_proteincoding_bed= chess_bed_file[chess_bed_file[12].isin(CHESS3_unique_gene_identifiers)].reset_index(drop=True)
	#subset to only include transcripts mapped to primary chromosomes
	unique_proteincoding_bed= unique_proteincoding_bed[-unique_proteincoding_bed[0].str.count("Un")==0].reset_index(drop=True)
	
	#Only keep primary transcript for each gene 
	#get transcript number within the gene
	unique_proteincoding_bed["transcript"]=unique_proteincoding_bed[3].apply(lambda x: int(x.split(".")[2]))
	#sort all CHESS transcripts by transcript number (within each gene)
	unique_proteincoding_bed=unique_proteincoding_bed.sort_values(by="transcript").reset_index(drop=True)
	#only keep first occuring transcript for each gene
	unique_proteincoding_bed=unique_proteincoding_bed.drop_duplicates(subset=12,keep="first").reset_index(drop=True)

	#Get exon ranges for each transcript
	#convert exon start and length positions to lists
	unique_proteincoding_bed["exon_start_pos"]=unique_proteincoding_bed[11].apply(lambda x: [int(i) for i in x.split(",") if len(i)>0])
	unique_proteincoding_bed["exon_length"]=unique_proteincoding_bed[10].apply(lambda x: [int(i) for i in x.split(",") if len(i)>0])
	#add starting position of transcript to get global coordinate start positions
	unique_proteincoding_bed["global_exon_start_pos"]=unique_proteincoding_bed.apply(lambda x: [i+x[1] for i in x.exon_start_pos],axis=1) 
	#add exon length to start position to get end position
	unique_proteincoding_bed["global_exon_end_pos"]=unique_proteincoding_bed.apply(lambda x: [start+length for start,length in zip(x.global_exon_start_pos,x.exon_length)],axis=1)
	#obtain start, end as range notation
	unique_proteincoding_bed["global_exon_ranges"]=unique_proteincoding_bed.apply(lambda x: [f"{start}-{end}" for start, end in zip(x.global_exon_start_pos, x.global_exon_end_pos)],axis=1)

	#Getting proper notation for CRISPick manual entry 
	#separate ranges with semicolon
	unique_proteincoding_bed['global_exon_ranges']=unique_proteincoding_bed['global_exon_ranges'].apply(lambda x: ";".join(x))
	#add in chromosome identifier
	unique_proteincoding_bed['global_exon_ranges']=unique_proteincoding_bed[0]+":-:"+unique_proteincoding_bed["global_exon_ranges"]
	#export file with just coordinates, identifiers
	unique_proteincoding_bed[['global_exon_ranges']].to_csv("Target-Gene-Identification-Data/CHESS3.1.3_unique_protein_coding_genes_coordinates.csv",header=False,index=False)

	#associate each exon range with CHESS identifier to associate guides back after CRISPick run
	unique_proteincoding_bed[[12,'global_exon_ranges']].to_csv("Target-Gene-Identification-Data/CHESS3.1.3_unique_protein_coding_genes_coordinates_with_identifier.csv",header=False,index=False)

if __name__ == "__main__":
    chess_unique_genes=get_unique_chess_genes()
    get_chess_coordinates(chess_unique_genes)

