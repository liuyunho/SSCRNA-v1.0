# SSCRNA-v1.0
single cell sequencing simulation program  

## Working Directory
Set up a working directory (WorkPath/), Create several new folders (fragment_data, pcr_control_data, umi_data, Simulation_Path).  
Put the files (Test_Data.txt, Error_Profile.txt, Quality_Profile.txt) in the working directory.  

## Running the main file (Main.py)
As recorded in Main.py file:  

#### Import packages:
import sys  
sys.path.append('WorkPath/')  
import SSCRNA as sc  
import time  
import os  
import pickle  

#### read the cell specific profile
filename = 'WorkPath/'+'/Test_Data.txt'  
#read the gene names and mRNA number  
cell_pro,gene_na=sc.Read_Gene_Table(filename)  

#### get the seq of genes (f_seq variable)
#sequences data can be found in ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.pc_transcripts.fa.gz  
f=r"gencode.v32.pc_transcripts.fa"  
#extract the gene names of transcripts fasta file  
f_na=sc.read_names_from_file(f)  
#split the name elements  
genes_in=sc.split_GENENAMES(f_na)  
#extract the Ensemble gene ids  
gene_name_list=sc.get_iterm_list(genes_in,2)  
#delete the ensemble version  
gene_name_list=sc.del_version(gene_name_list)  
sel_genes=sc.del_version(gene_na)  
#select genes corresponding ground truth (Test_Data.txt)  
sel_index=sc.get_gene_index(gene_name_list,sel_genes)  
#extract the corresponding sequences  
f_seq=sc.read_seqs_from_file(f,sel_index)  

#### get cell profile
f_na=sc.get_sel_gene_name(gene_name_list,sel_index)  
cell_pro_n=sc.get_map_gene_profile(cell_pro,sel_index)  

#### seperate cell name and profile to two list
cells_profile=[]  
cell_names=[]  
for each in cell_pro_n:  
	cell_names.append(each)  
	cells_profile.append(cell_pro_n[each])  

#### construct simulation sequencing library
#make real cell fragment database  
frag_list_path='WorkPath/'+'fragment_data/'  
sc.multi_cell2(frag_list_path,cells_profile,f_na,f_seq)  
#make pcr database  
pcr_list_path='WorkPath/'+'pcr_control_data/'  
sc.PCR_database(pcr_list_path,frag_list_path,3)  
#make UMI database  
umi_path='WorkPath/'+'umi_data/'  
sc.get_UMI_bank(frag_list_path,umi_path,8)  
#make barcode library:get_barcord_bank(length,least_distance,number)  
cell_barcord=sc.get_barcord_bank(27,2,22)  
#### save some data in pickle data
xc=[cell_barcord, f_na, f_seq]
with open('WorkPath/'+'cell_barcord_gene_name_seq.pickle', 'wb') as handle: 
    pickle.dump(xc, handle)

