#import packages
import sys
sys.path.append('WorkPath/')
import SSCRNA as sc
import time
import os
import pickle
#read the cell specific profile
filename = 'WorkPath/'+'/Test_Data.txt'
cell_pro,gene_na=sc.Read_Gene_Table(filename)
#get the sequences of the genes
f='WorkPath/'+"gencode.v32.pc_transcripts.fa"
f_na=sc.read_names_from_file(f)
genes_in=sc.split_GENENAMES(f_na)
gene_name_list=sc.get_iterm_list(genes_in,2)
gene_name_list=sc.del_version(gene_name_list)
sel_genes=sc.del_version(gene_na)
sel_index=sc.get_gene_index(gene_name_list,sel_genes)
f_seq=sc.read_seqs_from_file(f,sel_index)
#correct cell profile
f_na=sc.get_sel_gene_name(gene_name_list,sel_index)
cell_pro_n=sc.get_map_gene_profile(cell_pro,sel_index)
#define how many cell
cells_profile=[]
cell_names=[]
for each in cell_pro_n:
	cell_names.append(each)
	cells_profile.append(cell_pro_n[each])
#make real cell fragment database
frag_list_path='WorkPath/'+'fragment_data/'
sc.multi_cell2(frag_list_path,cells_profile,f_na,f_seq)
#make pcr database
pcr_list_path='WorkPath/'+'pcr_control_data/'
sc.PCR_database(pcr_list_path,frag_list_path,3)
real_size,PCR_list=sc.get_search_key(pcr_list_path)
#real_size:3412606400
#make UMI database
umi_path='WorkPath/'+'umi_data/'
sc.get_UMI_bank(frag_list_path,umi_path,8)
#make barcode library:get_barcord_bank(length,least_distance,number)
cell_barcord=sc.get_barcord_bank(27,2,22)
#common parameter#
xc=[cell_barcord, f_na, f_seq]
with open('WorkPath/'+'cell_barcord_gene_name_seq.pickle', 'wb') as handle: 
    pickle.dump(xc, handle)

