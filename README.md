# SSCRNA-v1.0
single cell sequencing simulation program  

## Working Directory
Set up a working directory (WorkPath/), Create several new folders (fragment_data, pcr_control_data, umi_data, Simulation_Path).  
Put the files (Test_Data.txt, Error_Profile.txt, Quality_Profile.txt) in the working directory.

## Running the main file (Main.py)
As recorded in Main.py file:

### Import packages:
import sys  
sys.path.append('WorkPath/')  
import SSCRNA as sc  
import time  
import os  
import pickle  

### read the cell specific profile
filename = r'/Test_Data_PATH/Test_Data.txt'
cell_pro,gene_na=sc.Read_Gene_Table(filename)

### get the seq of genes (f_seq variable)
### the refference seq download
#ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.pc_transcripts.fa.gz
f=r"gencode.v32.pc_transcripts.fa"
f_na=sc.read_names_from_file(f)
genes_in=sc.split_GENENAMES(f_na)
gene_name_list=sc.get_iterm_list(genes_in,2)
gene_name_list=sc.del_version(gene_name_list)
sel_genes=sc.del_version(gene_na)
sel_index=sc.get_gene_index(gene_name_list,sel_genes)
f_seq=sc.read_seqs_from_file(f,sel_index)

## correct cell profile
f_na=sc.get_sel_gene_name(gene_name_list,sel_index)
cell_pro_n=sc.get_map_gene_profile(cell_pro,sel_index)

## construct library
### make real cell fragment database
frag_list_path=r'/data3/users/liuyunhe/simulation_database_1022/real_break_data/'
sc.multi_cell2(frag_list_path,cells_profile,f_na,f_seq)
### make pcr database
pcr_list_path=r'/data3/users/liuyunhe/simulation_database_1022/pcr_control_data/'
sc.PCR_database(pcr_list_path,frag_list_path,3)
real_size,PCR_list=sc.get_search_key(pcr_list_path)
#real_size:3412606400
### make UMI database
umi_path=r'/data3/users/liuyunhe/simulation_database_1022/umi_data/'
sc.get_UMI_bank(frag_list_path,umi_path,8)

