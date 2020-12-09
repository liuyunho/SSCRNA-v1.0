import re as re
import random
import numpy as np
import linecache
import os
import time
from collections import Counter
import threading
import queue
def Read_Gene_Table(FILEPATH):
    cell_profiles={}
    gene_names=[]
    fr = open(FILEPATH, 'r')
    line=fr.readline()
    cell_a=line.replace('\n','').split('\t')
    for each in cell_a:
        cell_profiles.update({each: []})
    for line in fr:
        sub_a=line.replace('\n','').split('\t')
        for i in range(len(sub_a)):
            if i==0:
                gene_names.append(sub_a[i])
            else:
                cell_profiles[cell_a[i-1]].append(int(float((sub_a[i]))))
    return cell_profiles,gene_names
def read_names_from_file(FILENAME):
    fr = open(FILENAME, 'r')
    names=[] #define seq name list
    for line in fr:
        if line.startswith('>'): #decide if the start
            names.append(line.replace('>','').replace('\n',''))
    return names
def split_GENENAMES(GENE_NAMES):
    gene_infos=[]
    for each in GENE_NAMES:
        gene_infos.append(each.split('|'))
    return gene_infos
def del_version(GENE_LIST):
    gene_list=[]
    strinfo = re.compile(r'\..*')
    for each in GENE_LIST:
        gene_list.append(strinfo.sub('',each))
    return gene_list
def get_iterm_list(GENE_INFOS,INDEX):
    Iterm_List=[]
    for each in GENE_INFOS:
        Iterm_List.append(each[INDEX-1])
    return Iterm_List
def get_gene_index(GENE_LIST,SEL_GENES):
    Gene_Index=[]
    for each in SEL_GENES:
        if each in GENE_LIST:
            Gene_Index.append(GENE_LIST.index(each))
        else:
            Gene_Index.append('no')
    return Gene_Index
def read_seqs_from_file(FILENAME,SEL_INDEX):
    fr = open(FILENAME, 'r')
    index=(-1) #cause this part index start with 0
    seq_dic={}
    for line in fr:
        if line.startswith('>'): #decide if the start
            index=index+1
            if index in SEL_INDEX:
                seq_dic[index]=''
            continue #if the line start with >,no seq line,no need the following process
        if index in SEL_INDEX:
            seq_dic[index]+=line.replace('\n','')
    seqs = [] #define seq list
    for each in SEL_INDEX: #rank as the sel index order
        if each!='no':
            seqs.append(seq_dic[each])
    return seqs
def get_sel_gene_name(GENE_LIST,SEL_INDEX):
    names=[]
    for each in SEL_INDEX:
        if each!='no':
            names.append(GENE_LIST[each])
    return names
def get_map_gene_profile(CELL_PROFILES,SEL_INDEX):
    index=[]
    CELL_PROFILES_new={}
    for i in range(len(SEL_INDEX)):
        if SEL_INDEX[i]!='no':
            index.append(i)
    for each in CELL_PROFILES:
        CELL_PROFILES_new[each]=[CELL_PROFILES[each][each_i] for each_i in index]
    return CELL_PROFILES_new
def get_fragment_count(Real_Profile_in_a_Cell,Seq_Sequences):
    length_list=[] #a list contain the seqs' length
    for each in Seq_Sequences:
        length_list.append(len(each))
    length_ratio=[round(each/min(length_list)) for each in length_list] #get the length ratio,means the shortest length as 1
    fragment_count=list(map(lambda a:a[0]*a[1],zip(length_ratio,Real_Profile_in_a_Cell))) #the production of length ratio and real profile
    return fragment_count
def random_select(LENGTH,FRAG_SIZE):
    if LENGTH>FRAG_SIZE:
        return random.choice(range(1,(LENGTH-FRAG_SIZE+1)))
    else: #if length == frag size, range function return null list
        return 1 #if length is samller than length, then return the first position
def break_seq(LIBRARY_PATH,CELL_INDEX,Real_Profile_in_a_Cell,SEQ_NAMES,SEQ_SEQUENCES,Fragment_Size=250,Fragment_Dv=25):
    correct_fragment_count_in_a_cell=get_fragment_count(Real_Profile_in_a_Cell,SEQ_SEQUENCES) #correct fragment count by seqs' length
    fr=open(LIBRARY_PATH+str(CELL_INDEX),'w')
    for i in range(len(SEQ_NAMES)): #for each gene
        seq_leng=len(SEQ_SEQUENCES[i]) #each seq length
        sub_length=[round(random.normalvariate(Fragment_Size, Fragment_Dv)) for k in range(correct_fragment_count_in_a_cell[i])] #the fragment count is the real gene count, for every fragment, we random define the fragment size
        sub_length=[each if each<=seq_leng else seq_leng for each in sub_length] #fragment size cannot larger than gene size
        sub_start=[random_select(seq_leng,each) for each in sub_length] #for each fragment, we random choose a start position in the gene
        sub_end=list(np.array(sub_start)+np.array(sub_length)-1) #according to the fragment size dedined as above, we get the end point
        sub_fragment=list(zip(sub_start, sub_end)) #get every start and end to a tuple
        #write to file
        fr.writelines([SEQ_NAMES[i],'\t'])
        if len(sub_fragment)==0:
            fr.writelines(['\n'])
            continue
        for j in range(len(sub_fragment)):
            sub_a=[str(each) for each in list(sub_fragment[j])]
            sub_a.insert(1,',')
            fr.writelines(sub_a)
            if j==(len(sub_fragment)-1):
                fr.writelines(['\n'])
            else:
                fr.writelines([';'])
def break_seq2(LIBRARY_PATH,CELL_INDEX,Real_Profile_in_a_Cell,SEQ_NAMES,SEQ_SEQUENCES,Delete_Size=50,Delete_Dv=25):
    correct_fragment_count_in_a_cell=Real_Profile_in_a_Cell #dont need correct
    fr=open(LIBRARY_PATH+str(CELL_INDEX),'w')
    for i in range(len(SEQ_NAMES)): #for each gene
        seq_leng=len(SEQ_SEQUENCES[i]) #each seq length
        sub_length=[(seq_leng-round(random.normalvariate(Delete_Size, Delete_Dv))) for k in range(correct_fragment_count_in_a_cell[i])] #the fragment count is the real gene count, for every fragment, we random define the fragment size
        sub_length=[each if each<=seq_leng else seq_leng for each in sub_length] #fragment size cannot larger than gene size
        sub_start=[random_select(seq_leng,each) for each in sub_length] #for each fragment, we random choose a start position in the gene
        sub_end=list(np.array(sub_start)+np.array(sub_length)-1) #according to the fragment size dedined as above, we get the end point
        sub_fragment=list(zip(sub_start, sub_end)) #get every start and end to a tuple
        #write to file
        fr.writelines([SEQ_NAMES[i],'\t'])
        if len(sub_fragment)==0:
            fr.writelines(['\n'])
            continue
        for j in range(len(sub_fragment)):
            sub_a=[str(each) for each in list(sub_fragment[j])]
            sub_a.insert(1,',')
            fr.writelines(sub_a)
            if j==(len(sub_fragment)-1):
                fr.writelines(['\n'])
            else:
                fr.writelines([';'])
def multi_cell(LIBRARY_PATH,Cell_Profile_List,SEQ_NAMES,SEQ_SEQUENCES,Fragment_Size=250,Fragment_Dv=25):
    for i in range(len(Cell_Profile_List)): #for every cell fragment profile, we get the fragment list in a cell and combine into fragment list
        break_seq(LIBRARY_PATH,(i+1),Cell_Profile_List[i],SEQ_NAMES,SEQ_SEQUENCES,Fragment_Size,Fragment_Dv)
def multi_cell2(LIBRARY_PATH,Cell_Profile_List,SEQ_NAMES,SEQ_SEQUENCES,Delete_Size=50,Delete_Dv=25):
    for i in range(len(Cell_Profile_List)): #for every cell fragment profile, we get the fragment list in a cell and combine into fragment list
        break_seq2(LIBRARY_PATH,(i+1),Cell_Profile_List[i],SEQ_NAMES,SEQ_SEQUENCES,Delete_Size,Delete_Dv)
def get_line_context(file_path, line_number):
    return linecache.getline(file_path, line_number).strip()
def split_gene_fragment(GENE_Fragment):
    GENE_Fragment=GENE_Fragment.split('\t')
    if GENE_Fragment[1]=='':
        return GENE_Fragment[0],[]
    else:
        fragments=GENE_Fragment[1].split(';')
        return GENE_Fragment[0],[[int(eeach) for eeach in each.split(',')] for each in fragments]
def show_fragment(LIBRARY_PATH,SEQ_SEQUENCES,Cell_Index=1,Gene_Index=1,Fragment_Index=1):
    Fragment_Index=Fragment_Index-1 #change the index to python index
    seq_name,fragement_list=split_gene_fragment(get_line_context(LIBRARY_PATH+str(Cell_Index),Gene_Index))
    fragment_seq=SEQ_SEQUENCES[Gene_Index-1][fragement_list[Fragment_Index][0]:(1+fragement_list[Fragment_Index][1])]
    print(seq_name) #print seq name
    print(fragment_seq) #print seq
    print(len(fragment_seq)) #print seq length
def list_insert(LIST,SEP):
    LIST_n=[]
    for i in range(len(LIST)):
        LIST_n.append(LIST[i])
        LIST_n.append(SEP)
    LIST_n.pop()
    return LIST_n
def PCR_database(PCR_CONTROL_PATH,LIBRARY_PATH,CIRCLE=2):
    cell_list=os.listdir(LIBRARY_PATH)
    for each in cell_list:
        fr=open(LIBRARY_PATH+each,'r')
        fw=open(PCR_CONTROL_PATH+each,'w')
        for line in fr:
            line_na,line_frag=split_gene_fragment(line.replace('\n',''))
            fw.writelines([line_na,'\t'])
            if len(line_frag)==0:
                fw.writelines(['\n']) #if the real count is 0,then the space is empty
                continue
            pcr=[str(2**CIRCLE) for eeach in line_frag]
            pcr=list_insert(pcr,',')
            pcr.extend(['\n'])
            fw.writelines(pcr)
        fr.close()
        fw.close()
def edit_pcr_control(PCR_CONTROL_PATH,NEW_VALUE,Cell_Index=1,Gene_Index=1,Fragment_Index=1):
    fr=open(PCR_CONTROL_PATH+str(Cell_Index),'r')
    lines=fr.readlines()
    fr.close()
    sub_line=lines[Gene_Index-1]
    sub_line=sub_line.replace('\n','').split('\t')
    sub_frag=sub_line[1].split(',')
    sub_frag[Fragment_Index-1]=str(NEW_VALUE)
    sub_line=[sub_line[0],'\t']
    sub_line.extend(list_insert(sub_frag,','))
    sub_line.append('\n')
    lines[Gene_Index-1]=''.join(sub_line)
    fw=open(PCR_CONTROL_PATH+str(Cell_Index),'w')
    fw.writelines(lines)
    fw.close()
def get_count_sum(PCR_Fragment):
    pcr_list=PCR_Fragment.replace('\n','').split('\t')
    if pcr_list[1]=='':  #if no fragment
        return 0
    else:
        sub_pcr_frag=[int(each) for each in pcr_list[1].split(',')]
        return sum(sub_pcr_frag)
def get_search_key(PCR_CONTROL_PATH):
    cell_list=[int(each) for each in os.listdir(PCR_CONTROL_PATH)]
    cell_list.sort()
    pcr_list=[]
    database_size=0
    for i in cell_list:
        sub_pcr_list=[]
        fr=open(PCR_CONTROL_PATH+str(i),'r')
        for line in fr:
            sub_sum=get_count_sum(line)
            database_size=database_size+sub_sum
            sub_pcr_list.append(sub_sum)
        fr.close()
        pcr_list.append(sub_pcr_list)
    return database_size,pcr_list
def get_cell_count(PCR_LIST):
    return [sum(each) for each in PCR_LIST]
def find_expend_index(Database_Expend,Random_Number):
    Database_Expend=np.array(Database_Expend)
    index_i=1
    index_j=len(Database_Expend)
    while (index_j-index_i)>1: #the problem is: i cannot be higher(can be a loop:j=i), or j cannot be lower(can be a loop:j-i=1)
        index=int((index_i+index_j)/2) #the middle one(.5 to .0:only retain the integer)
        if np.sum(Database_Expend[:index])<Random_Number:
            index_i=index
        elif np.sum(Database_Expend[:index])>=Random_Number:
            index_j=index
    if np.sum(Database_Expend[:index_i])>=Random_Number:
        return index_i
    else:
        return index_j
def get_frag_index(PCR_CONTROL_PATH,PCR_LIST,Random_Number):
    #find cell index
    cell_count_sum=get_cell_count(PCR_LIST)
    cell_index=find_expend_index(cell_count_sum,Random_Number)
    #find gene index
    Random_Number=Random_Number-sum(cell_count_sum[:(cell_index-1)])
    gene_count_sum=PCR_LIST[cell_index-1] #cause the python list index rule
    gene_index=find_expend_index(gene_count_sum,Random_Number)
    #find fragment index
    Random_Number=Random_Number-sum(gene_count_sum[:(gene_index-1)])
    fragment_line=get_line_context(PCR_CONTROL_PATH+str(cell_index),gene_index)
    fragment_count_sum=[int(each) for each in fragment_line.split('\t')[1].split(',')] #must not be empty fragment count
    fragment_index=find_expend_index(fragment_count_sum,Random_Number)
    return [cell_index,gene_index,fragment_index]
def cultimate_list(LIST):
    CULM_LIST=[]
    for i in range(len(LIST)):
        if i==0:
            CULM_LIST.append(LIST[i])
            continue
        CULM_LIST.append(sum([CULM_LIST[(i-1)],LIST[i]]))
    return CULM_LIST
def find_cumulate_index(Cumulative_List,Sorted_Index_List):
    Index=[]
    index_j=0
    index_i=0
    while index_j<len(Cumulative_List):
        while index_i<len(Sorted_Index_List):
            if Sorted_Index_List[index_i]<=Cumulative_List[index_j]:
                Index.append(index_j+1)
                index_i=index_i+1
            else:
                break
        index_j=index_j+1
    return Index
def find_cell_index(Cumulative_List,Sorted_Index_List):
    Index=[]
    index_j=0
    index_i=0
    while index_j<len(Cumulative_List):
        while index_i<len(Sorted_Index_List):
            if Sorted_Index_List[index_i]<=Cumulative_List[index_j]:
                Index.append(index_j+1)
                index_i=index_i+1
            else:
                break
        print('\r'+'*'*(((index_j*100)//(len(Cumulative_List)-1))//2)+str(((index_j*100)//(len(Cumulative_List)-1)))+'%', end='')
        index_j=index_j+1
    return Index
def find_gene_index(PCR_Cell_LIST,Cell_Count_cumuLIST,Cell_Index_LIST,LIBRARY_INDEX):
    sel_cell=list(Counter(Cell_Index_LIST).keys())
    sel_cell.sort()
    ##############
    gene_index=[]
    index=0
    for each_cell in sel_cell:
        gene_count_culm=cultimate_list(PCR_Cell_LIST[each_cell-1])
        #get the sub library index, each sub index means of in one cell
        if each_cell==sel_cell[-1]:
            if each_cell==1:
                sub_lib_index=LIBRARY_INDEX[slice(Cell_Index_LIST.index(sel_cell[index]),len(Cell_Index_LIST))]
            else:
                sub_lib_index=[(ii-Cell_Count_cumuLIST[each_cell-2]) for ii in LIBRARY_INDEX[slice(Cell_Index_LIST.index(sel_cell[index]),len(Cell_Index_LIST))]]
        elif each_cell==sel_cell[-1]: #equal to the last cell in cell index list, need to first decide
            sub_lib_index=[(ii-Cell_Count_cumuLIST[each_cell-2]) for ii in LIBRARY_INDEX[slice(Cell_Index_LIST.index(sel_cell[index]),len(Cell_Index_LIST))]]
        elif each_cell==1:
            sub_lib_index=LIBRARY_INDEX[slice(Cell_Index_LIST.index(sel_cell[index]),Cell_Index_LIST.index(sel_cell[index+1]))]
        else:
            sub_lib_index=[(ii-Cell_Count_cumuLIST[each_cell-2]) for ii in LIBRARY_INDEX[slice(Cell_Index_LIST.index(sel_cell[index]),Cell_Index_LIST.index(sel_cell[index+1]))]]
        gene_index.extend(find_cumulate_index(gene_count_culm,sub_lib_index))
        index=index+1
        #print the progress
        print('\r'+'*'*(((each_cell*100)//(sel_cell[-1]))//2)+str(((each_cell*100)//(sel_cell[-1])))+'%', end='')
    return gene_index
def find_fragment_index_in_cell(PCR_LIST,PCR_CONTROL_PATH,CELL,SUB_LIBRARY_INDEX,SUB_GENE_INDEX):
    gene_index=[]
    sel_gene=[]
    with open(PCR_CONTROL_PATH+str(CELL),'r') as fr:
        fr_index=1
        for each in fr:
            if fr_index in SUB_GENE_INDEX:
                gene_index.append([int(ii) for ii in each.split('\t')[1].replace('\n','').split(',')])
                sel_gene.append(fr_index)
            fr_index=fr_index+1
    gene_count_culm=cultimate_list(PCR_LIST[CELL-1])
    fragment_index=[]
    #SUB_LIBRARY_index=lib_index[slice(b.index(1),b.index(2))]
    #SUB_GENE_index=aaa[slice(b.index(1),b.index(2))]
    index=0
    for i in sel_gene:
        if i==sel_gene[-1]:
            if i==1:
                sub_lib_index=SUB_LIBRARY_INDEX[slice(SUB_GENE_INDEX.index(sel_gene[index]),(len(SUB_GENE_INDEX)))]
            else:
                sub_lib_index=[ii-gene_count_culm[i-2] for ii in SUB_LIBRARY_INDEX[slice(SUB_GENE_INDEX.index(sel_gene[index]),(len(SUB_GENE_INDEX)))]]
        elif i==1:
            sub_lib_index=SUB_LIBRARY_INDEX[slice(SUB_GENE_INDEX.index(sel_gene[index]),SUB_GENE_INDEX.index(sel_gene[index+1]))]
        elif i==sel_gene[-1]:
            sub_lib_index=[ii-gene_count_culm[i-2] for ii in SUB_LIBRARY_INDEX[slice(SUB_GENE_INDEX.index(sel_gene[index]),(len(SUB_GENE_INDEX)))]]
        else:
            sub_lib_index=[ii-gene_count_culm[i-2] for ii in SUB_LIBRARY_INDEX[slice(SUB_GENE_INDEX.index(sel_gene[index]),SUB_GENE_INDEX.index(sel_gene[index+1]))]]
        fragment_index.extend(find_cumulate_index(cultimate_list(gene_index[index]),sub_lib_index))
        index=index+1
    return fragment_index
def find_fragment_index(PCR_LIST,PCR_CONTROL_PATH,LIBRARY_INDEX,CELL_INDEX,GENE_INDEX,cell_count_culm):
    sel_cell=[]
    for each in CELL_INDEX:
        if each not in sel_cell:
            sel_cell.append(each)
    all_frag_index=[]
    index=0
    for each_cell in sel_cell:
        if each_cell==sel_cell[-1]:
            if each_cell==1:
                SUB_LIBRARY_index=LIBRARY_INDEX[slice(CELL_INDEX.index(sel_cell[index]),len(CELL_INDEX))]
                SUB_GENE_index=GENE_INDEX[slice(CELL_INDEX.index(sel_cell[index]),len(CELL_INDEX))]
            else:
                SUB_LIBRARY_index=[(ii-cell_count_culm[each_cell-2]) for ii in LIBRARY_INDEX[slice(CELL_INDEX.index(sel_cell[index]),len(CELL_INDEX))]]
                SUB_GENE_index=GENE_INDEX[slice(CELL_INDEX.index(sel_cell[index]),len(CELL_INDEX))]
        elif each_cell==sel_cell[-1]:
            SUB_LIBRARY_index=[(ii-cell_count_culm[each_cell-2]) for ii in LIBRARY_INDEX[slice(CELL_INDEX.index(sel_cell[index]),len(CELL_INDEX))]]
            SUB_GENE_index=GENE_INDEX[slice(CELL_INDEX.index(sel_cell[index]),len(CELL_INDEX))]
        elif each_cell==1:
            SUB_LIBRARY_index=LIBRARY_INDEX[slice(CELL_INDEX.index(sel_cell[index]),CELL_INDEX.index(sel_cell[index+1]))]
            SUB_GENE_index=GENE_INDEX[slice(CELL_INDEX.index(sel_cell[index]),CELL_INDEX.index(sel_cell[index+1]))]
        else:
            SUB_LIBRARY_index=[(ii-cell_count_culm[each_cell-2]) for ii in LIBRARY_INDEX[slice(CELL_INDEX.index(sel_cell[index]),CELL_INDEX.index(sel_cell[index+1]))]]
            SUB_GENE_index=GENE_INDEX[slice(CELL_INDEX.index(sel_cell[index]),CELL_INDEX.index(sel_cell[index+1]))]
        all_frag_index.extend(find_fragment_index_in_cell(PCR_LIST,PCR_CONTROL_PATH,each_cell,SUB_LIBRARY_index,SUB_GENE_index))
        ##############
        index=index+1
        print('\r'+'*'*(((each_cell*100)//sel_cell[-1])//2)+str(((each_cell*100)//sel_cell[-1]))+'%', end='')
    return all_frag_index
def get_all_Count_UMI(PCR_CONTROL_PATH,PCR_LIST,LIBRARY_INDEX):
    #
    print('Getting Cells cumulative count')
    cell_count_sum=get_cell_count(PCR_LIST)
    cell_count_culm=cultimate_list(cell_count_sum)
    print('Get '+str(len(cell_count_culm))+' cells')
    print('The real library size is:  '+str(cell_count_culm[len(cell_count_culm)-1]))
    print('------Have Done!------')
    #
    print('\nSorting the random index library')
    LIBRARY_INDEX.sort()
    print('The range of the index library is: ['+str(min(LIBRARY_INDEX))+', ',str(max(LIBRARY_INDEX))+']')
    print('------Have Done!------')
    #cell index
    print('\nProcessing cell index calculation')
    cell_index_list=find_cell_index(cell_count_culm,LIBRARY_INDEX)
    print('\nThe range of the index library is: ['+str(min(cell_index_list))+', ',str(max(cell_index_list))+']')
    print('------Have Done!------')
    #gene index
    print('\nProcessing gene index calculation')
    gene_index_list=find_gene_index(PCR_LIST,cell_count_culm,cell_index_list,LIBRARY_INDEX)
    print('\n------Have Done!------')
    #fragment index
    print('\nProcessing fragment index calculation')
    fragment_index_list=find_fragment_index(PCR_LIST,PCR_CONTROL_PATH,LIBRARY_INDEX,cell_index_list,gene_index_list,cell_count_culm)
    print('\n------Have Done!------')
    #
    return [list(each) for each in zip(cell_index_list,gene_index_list,fragment_index_list)]
def random_seq(LENGTH):
    return ''.join([random.choice(['A','T','G','C']) for i in range(LENGTH)])
def seq_distance(SEQ_1,SEQ_2):
    index=0
    for i in range(len(SEQ_1)):
        if SEQ_1[i]!=SEQ_2[i]:
            index=index+1
    return index
def get_new_barcord(LENGTH,MIN_DISTANT,Barcord_Bank):
    random_barcord=random_seq(LENGTH)
    dis_in_seqs=seq_distance(random_barcord,Barcord_Bank[0])
    if len(Barcord_Bank)==1:
        if dis_in_seqs>=MIN_DISTANT:
            return random_barcord
        else:
            return get_new_barcord(LENGTH,MIN_DISTANT,Barcord_Bank) #if not size up, then get a new one
    else:
        i=1
        while (dis_in_seqs>=MIN_DISTANT)&(i<len(Barcord_Bank)):
            dis_in_seqs=min([dis_in_seqs,seq_distance(random_barcord,Barcord_Bank[i])])
            i=i+1
        if dis_in_seqs>=MIN_DISTANT:
            return random_barcord
        else:
            return get_new_barcord(LENGTH,MIN_DISTANT,Barcord_Bank) #if not size up with the barcord in the bank, then get a new one
def get_barcord_bank(LENGTH,MIN_DISTANT,BANK_SIZE):
    barcord_bank=[random_seq(LENGTH)] #get the first one
    i=1
    while i<BANK_SIZE:
        barcord_bank.append(get_new_barcord(LENGTH,MIN_DISTANT,barcord_bank))
        i=i+1
    return barcord_bank
def get_UMI_bank(LIBRARY_PATH,UMI_PATH,UMI_Length):
    cells=os.listdir(LIBRARY_PATH)
    for each in cells:
        fr=open(LIBRARY_PATH+each,'r')
        fw=open(UMI_PATH+each,'w')
        for line in fr:
            line_na,line_frag=split_gene_fragment(line.replace('\n',''))
            fw.writelines([line_na,'\t'])
            if len(line_frag)==0:
                fw.writelines(['\n']) #if the real count is 0,then the space is empty
                continue
            umi_list=[''.join([random.choice(['0','1','2','3']) for i in range(UMI_Length)]) for eeach in line_frag]
            umi_list=list_insert(umi_list,',')
            umi_list.extend(['\n'])
            fw.writelines(umi_list)
        fr.close()
        fw.close()
def get_UMI_barcord(UMI_PATH,Full_Index):
    umi_line=get_line_context(UMI_PATH+str(Full_Index[0]),Full_Index[1])
    umi_line=umi_line.split('\t')[1].split(',')
    return ''.join([['A','T','G','C'][int(b)] for b in list(umi_line[Full_Index[2]-1])])
def get_library_index(Database_Size,Library_Size):
    library_index=random.sample(range(Database_Size),Library_Size)
    return [i+1 for i in library_index]
def get_fragment_seq(LIBRARY_PATH,SEQ_SEQUENCES,Full_Index):
    fragment_index=Full_Index[2]-1 #change the index to python index
    Cell_Index=Full_Index[0]
    Gene_Index=Full_Index[1]
    seq_name,fragement_list=split_gene_fragment(get_line_context(LIBRARY_PATH+str(Cell_Index),Gene_Index))
    fragment_seq=SEQ_SEQUENCES[Gene_Index-1][fragement_list[fragment_index][0]:(1+fragement_list[fragment_index][1])]
    return fragment_seq
def get_complimentary_seq(SEQ_1):
    SEQ_2=[['T','A','C','G'][['A','T','G','C'].index(each)] for each in SEQ_1]
    SEQ_2=''.join(SEQ_2)
    return SEQ_2[::-1]
def get_full_seq(LIBRARY_PATH,SEQ_SEQUENCES,Full_Index,Cell_Barcord,UMI_PATH,ADAPTOR):
    cb=Cell_Barcord[Full_Index[0]-1]
    ub=get_UMI_barcord(UMI_PATH,Full_Index)
    return ADAPTOR+cb[:9]+'ACTGGCCTGCGA'+cb[9:18]+'GGTAGCGGCGACA'+cb[18:27]+ub+'T'+get_fragment_seq(LIBRARY_PATH,SEQ_SEQUENCES,Full_Index)+get_complimentary_seq(ADAPTOR)
def decide_where(Sparse_Count_List,LIBRARY_INDEX): #LIBRARY_INDEX:[cell index, gene index]
    if len(Sparse_Count_List)==0:
        return 'new'
    index_bank=[each[:2] for each in Sparse_Count_List]
    if LIBRARY_INDEX in index_bank:
        return index_bank.index(LIBRARY_INDEX)
    else:
        return 'new'
def get_one_read(FULL_SEQ,STRAND,START_RANGE=[0,10],LENGTH=125):
    start_position=random.choice(range((START_RANGE[0]),(START_RANGE[1])))
    end_position=min([(start_position+LENGTH),len(FULL_SEQ)])
    if STRAND==1:
        return FULL_SEQ[start_position:end_position]
    elif STRAND==2:
        return get_complimentary_seq(FULL_SEQ)[start_position:end_position]
    else:
        print('input the right strand')
def random_error():
    return random.randrange(0, 10, 1)/1000
def what_mutation(OLD_BASE,Random_Number,Mutant_Profile):
    Random_Number=Random_Number-(1-np.sum(np.array(Mutant_Profile))) #minus the probability no mutant
    index=0
    while Random_Number>0:
        Random_Number=Random_Number-Mutant_Profile[index] #minus the probability formar mutant
        index=index+1
    index=index-1
    if index<5: #delete or substitute mutant
        return ['','A','T','G','C'][index]
    else: #insert mutant
        return OLD_BASE+['','-','-','-','-','A','T','G','C'][index]
def get_base(OLD_BASE,POSITION,Error_Profile,STRAND):
    STRAND=STRAND-1
    POSITION=POSITION-1
    j=['A','T','G','C'].index(OLD_BASE)
    random_number=random.random() #random get a number 0-1
    if random_number<=(1-np.sum(np.array(Error_Profile[STRAND][j][POSITION]))):
        mutant_index=0
        return OLD_BASE,mutant_index
    else:
        mutant_index=1
        return what_mutation(OLD_BASE,random_number,Error_Profile[STRAND][j][POSITION]),mutant_index
def base_quality(Quality):
    Quality=Quality+round(random.normalvariate(0,0.5)) #a normal distribution with mean=0,sd=2
    if Quality<33:
        return chr(33)
    elif Quality>74:
        return chr(74)
    else:
        return chr(Quality)
def get_base_quality(NEW_BASE,POSITION,STRAND,Mutant_Index,Quality_Profile):
    if Mutant_Index==1:
        if NEW_BASE=='':
            return ''
        elif len(NEW_BASE)==2:
            return chr(33)+chr(33)
        else:
            return chr(33)
    else:
        STRAND=STRAND-1
        POSITION=POSITION-1
        return base_quality(Quality_Profile[STRAND][['A','T','G','C'].index(NEW_BASE)][POSITION])
def random_quality_profile(): #low-high-middle
    return [random.choice(range(33,51)) for i in range(25)]+[75 for i in range(75)]+[random.choice(range(40,66)) for i in range(25)]
def get_seq_read(OLD_SEQUENCE,STRAND,Error_Profile,Quality_Profile):
    new_sequence=''
    quality_score=''
    for i in range(len(OLD_SEQUENCE)):
        position=i+1
        new_sub,mutant_index=get_base(OLD_SEQUENCE[i],position,Error_Profile,STRAND)
        new_sequence=new_sequence+new_sub
        quality_score=quality_score+get_base_quality(new_sub,position,STRAND,mutant_index,Quality_Profile)
    return [new_sequence,quality_score]
def get_seq_title(SEQ_NAMES,Cell_Barcord,Full_Index,STRAND,Random_Index):
    if STRAND==1:
        return '@RDF'+str(Random_Index)+'_'+SEQ_NAMES[Full_Index[1]-1]+':'+Cell_Barcord[Full_Index[0]-1]
    else:
        return '@RDR'+str(Random_Index)+'_'+SEQ_NAMES[Full_Index[1]-1]+':'+Cell_Barcord[Full_Index[0]-1]
class myThread (threading.Thread):
    def __init__(self, name, q):
        threading.Thread.__init__(self,name=name)
        self.q = q
    def run(self):
        print (self.name+' start working now!')
        process_data(self.name, self.q) #in quene, will not stop
        print ('\n'+self.name+' have finish his work!')
def process_data(threadName,q):
    FILE=simulation_frag_to_file.FILE_PATH+'/'+simulation_frag_to_file.FILE_NAME+'_line'+threadName.split('-')[-1]
    flush_empty_index=0
    write_flush_1=[]
    write_flush_2=[]
    while not exitFlag:
        if not workQueue.empty():
            fragment_bed_full = q.get()
            full_seq=get_full_seq2(simulation_frag_to_file.SEQUENCES,simulation_frag_to_file.ADAPTOR,simulation_frag_to_file.Cell_Barcord,fragment_bed_full)
            read1_seq=get_one_read(full_seq,1,[0,0],125)
            read1_all=get_seq_read(read1_seq,1,simulation_frag_to_file.ERROR_Profile,simulation_frag_to_file.QUALITY_Profile)
            read2_seq=get_one_read(full_seq,2,[0,0],125)
            read2_all=get_seq_read(read2_seq,2,simulation_frag_to_file.ERROR_Profile,simulation_frag_to_file.QUALITY_Profile)
            #insert to flush
            write_flush_1.extend([get_seq_title2(simulation_frag_to_file.SEQUENCE_NAMES,simulation_frag_to_file.Cell_Barcord,fragment_bed_full,1),'\n'])
            write_flush_1.extend([read1_all[0],'\n'])
            write_flush_1.extend(['+','\n'])
            write_flush_1.extend([read1_all[1],'\n'])
            write_flush_2.extend([get_seq_title2(simulation_frag_to_file.SEQUENCE_NAMES,simulation_frag_to_file.Cell_Barcord,fragment_bed_full,2),'\n'])
            write_flush_2.extend([read2_all[0],'\n'])
            write_flush_2.extend(['+','\n'])
            write_flush_2.extend([read2_all[1],'\n'])
            flush_empty_index=1
            if len(write_flush_1)<1000:
                pass
            else:
                #queue up
                #queueLock.acquire()
                with open(FILE+'_R1.fq','a') as f1,open(FILE+'_R2.fq','a') as f2:
                    f1.writelines(write_flush_1)
                    f2.writelines(write_flush_2)
                #queueLock.release()
                #release off
                write_flush_1=[]
                write_flush_2=[]
                flush_empty_index=0
        elif flush_empty_index==1:
            #queueLock.acquire()
            with open(FILE+'_R1.fq','a') as f1,open(FILE+'_R2.fq','a') as f2:
                f1.writelines(write_flush_1)
                f2.writelines(write_flush_2)
            #queueLock.release()
            write_flush_1=[]
            write_flush_2=[]
            flush_empty_index=0
def simulation_frag_to_file(FILE_PATH,FILE_NAME,SEQUENCES,SEQUENCE_NAMES,ADAPTOR,Cell_Barcord,ERROR_Profile,QUALITY_Profile,LIBRARY_INDEX_Iterm,THREAD_Number):
    #share the variable
    simulation_frag_to_file.FILE_PATH=FILE_PATH
    simulation_frag_to_file.FILE_NAME=FILE_NAME
    simulation_frag_to_file.SEQUENCES=SEQUENCES
    simulation_frag_to_file.SEQUENCE_NAMES=SEQUENCE_NAMES
    simulation_frag_to_file.ADAPTOR=ADAPTOR
    simulation_frag_to_file.Cell_Barcord=Cell_Barcord
    simulation_frag_to_file.ERROR_Profile=ERROR_Profile
    simulation_frag_to_file.QUALITY_Profile=QUALITY_Profile
    simulation_frag_to_file.LIBRARY_INDEX_Iterm=LIBRARY_INDEX_Iterm
    #define global variable
    global exitFlag
    global workQueue
    global queueLock
    ################
    exitFlag= 0 
    queueLock = threading.Lock()
    threads = []
    #make a queue no limit
    workQueue = queue.Queue(-1)
    #make thread
    for tNum in range(THREAD_Number):
        thread = myThread(('Lisa Robort-'+str(tNum)), workQueue)
        thread.start()
        threads.append(thread)
    #stop the progress
    queueLock.acquire()
    #fill the queue
    for each_bed in LIBRARY_INDEX_Iterm:
        workQueue.put(each_bed)
    #get the full size
    size_queue=workQueue.qsize()
    #start release the queue
    queueLock.release()
    #check the progress every 2 seconds
    while not workQueue.empty():
        time.sleep(2)
        now_size=workQueue.qsize()
        print('\r'+'*'*((((size_queue-now_size)*100)//(size_queue))//2)+str((((size_queue-now_size)*100)//(size_queue)))+'%', end='')
    #pass the formar loop, stop the thread
    exitFlag = 1
    #waiting the full progress done
    for t in threads:
        t.join()
    print ("The simulation progress is done!")
def get_full_seq2(SEQ_SEQUENCES,ADAPTOR,Cell_Barcord,Fragment_bed_full):
    cell_index=Fragment_bed_full[0][0]
    gene_index=Fragment_bed_full[0][1]
    fragment_bed=Fragment_bed_full[1]
    cb=Cell_Barcord[cell_index-1]
    ub=''.join([['A','T','G','C'][int(b)] for b in list(Fragment_bed_full[2])])
    return cb+ub+'TTTTTTTTTTTTTTTTTT'+get_complimentary_seq(SEQ_SEQUENCES[gene_index-1][fragment_bed[0]:(1+fragment_bed[1])])
def get_seq_title2(SEQ_NAMES,Cell_Barcord,Fragment_bed_full,STRAND):
    if STRAND==1:
        return '@RD_'+str(Fragment_bed_full[1][0])+';'+str(Fragment_bed_full[1][1])+'_'+SEQ_NAMES[Fragment_bed_full[0][1]-1]+':'+Cell_Barcord[Fragment_bed_full[0][0]-1]
    else:
        return '@RD_'+str(Fragment_bed_full[1][0])+';'+str(Fragment_bed_full[1][1])+'_'+SEQ_NAMES[Fragment_bed_full[0][1]-1]+':'+Cell_Barcord[Fragment_bed_full[0][0]-1]
def split_gene_fragment2(GENE_Fragment,fragment_list):
    GENE_Fragment=GENE_Fragment.split('\t')
    fragments=[[int(eeach) for eeach in each.split(',')] for each in GENE_Fragment[1].split(';')]
    return [fragments[i-1] for i in fragment_list]
def split_UMI_fragment(UMI_line,fragment_list):
    UMI_line=UMI_line.replace('\n','').split('\t')
    fragments=UMI_line[1].split(',')
    return [fragments[i-1] for i in fragment_list]
def get_fragment_seq2(LIBRARY_PATH,UMI_PATH,Sub_Count_Sparse):
    #split the count sparse
    cell_index=[]
    gene_index=[]
    fragment_index=[]
    sel_cell=[]
    libary_size=len(Sub_Count_Sparse)
    print('Split the FULL INDEX:')
    for each in Sub_Count_Sparse:
        cell_index.append(each[0])
        gene_index.append(each[1])
        fragment_index.append(each[2])
        if each[0] not in sel_cell:
            sel_cell.append(each[0])
        #print('\r'+'▇'*(((len(cell_index)*100)//libary_size)//2)+str(((len(cell_index)*100)//libary_size))+'%', end='')
    print('\nHave Done!\n')
    fragment_bed=[]
    ##################################
    print('\nStart fragment bed searching:')
    for i in range(len(sel_cell)):
        #get sub gene index and fragment index in one cell
        if sel_cell[i]==sel_cell[-1]: #the last term. 
            sub_gene_index=gene_index[slice(cell_index.index(sel_cell[i]),len(cell_index))]
            sub_fragment_index=fragment_index[slice(cell_index.index(sel_cell[i]),len(cell_index))]
        else:
            sub_gene_index=gene_index[slice(cell_index.index(sel_cell[i]),cell_index.index(sel_cell[i+1]))]
            sub_fragment_index=fragment_index[slice(cell_index.index(sel_cell[i]),cell_index.index(sel_cell[i+1]))]
        #######################
        sel_gene=list(Counter(sub_gene_index).keys())
        sel_gene.sort()
        index_i=1 #indicate the line number and also the gene number
        with open(LIBRARY_PATH+str(sel_cell[i]),'r') as fr, open(UMI_PATH+str(sel_cell[i]),'r') as fu:
            index_j=0 #indict the sel_gene entry
            for line,line_u in zip(fr,fu):
                if index_i in sel_gene:
                    if index_i==sel_gene[-1]: #dont need break,cause the length of sel_cell is one
                        ssub_fragment_index=sub_fragment_index[slice(sub_gene_index.index(sel_gene[index_j]),len(sub_gene_index))]
                    else:
                        ssub_fragment_index=sub_fragment_index[slice(sub_gene_index.index(sel_gene[index_j]),sub_gene_index.index(sel_gene[index_j+1]))]
                    #fragment_bed.extend(split_gene_fragment2(line,ssub_fragment_index))
                    fragment_bed.extend(list(zip([[sel_cell[i],index_i] for n in range(len(ssub_fragment_index))],split_gene_fragment2(line,ssub_fragment_index),split_UMI_fragment(line_u,ssub_fragment_index))))
                    index_j=index_j+1
                index_i=index_i+1
        print('\r'+'*'*(((i*100)//(len(sel_cell)-1))//2)+str(((i*100)//(len(sel_cell)-1)))+'%', end='')
    print('\nHave Done!')
    return fragment_bed
