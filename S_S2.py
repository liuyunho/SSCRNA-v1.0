import re as re
import random
import numpy as np
import linecache
import os
import time
from collections import Counter
def get_full_seq2(SEQ_SEQUENCES,ADAPTOR,Cell_Barcord,Fragment_bed_full):
    cell_index=Fragment_bed_full[0][0]
    gene_index=Fragment_bed_full[0][1]
    fragment_bed=Fragment_bed_full[1]
    cb=Cell_Barcord[cell_index-1]
    ub=''.join([['A','T','G','C'][int(b)] for b in list(Fragment_bed_full[2])])
    return cb+ub+'TTTTTTTTTTTTTTTTTT'+get_complimentary_seq(SEQ_SEQUENCES[gene_index-1][fragment_bed[0]:(1+fragment_bed[1])])
def get_complimentary_seq(SEQ_1):
    SEQ_2=[['T','A','C','G'][['A','T','G','C'].index(each)] for each in SEQ_1]
    SEQ_2=''.join(SEQ_2)
    return SEQ_2[::-1]
#add 19 former is wrong
def get_seq_title2(SEQ_NAMES,Cell_Barcord,Fragment_bed_full,STRAND):
    if STRAND==1:
        return '@RD_'+str(Fragment_bed_full[1][0])+';'+str(Fragment_bed_full[1][1])+'_'+SEQ_NAMES[Fragment_bed_full[0][1]-1]+':'+Cell_Barcord[Fragment_bed_full[0][0]-1]
    else:
        return '@RD_'+str(Fragment_bed_full[1][0])+';'+str(Fragment_bed_full[1][1])+'_'+SEQ_NAMES[Fragment_bed_full[0][1]-1]+':'+Cell_Barcord[Fragment_bed_full[0][0]-1]
def get_one_read(FULL_SEQ,STRAND,START_RANGE=[0,1],LENGTH=125):
    start_position=random.choice(range((START_RANGE[0]),(START_RANGE[1])))
    end_position=min([(start_position+LENGTH),len(FULL_SEQ)])
    if STRAND==1:
        return FULL_SEQ[start_position:end_position]
    elif STRAND==2:
        return get_complimentary_seq(FULL_SEQ)[start_position:end_position]
    else:
        print('input the right strand')
def get_seq_read(OLD_SEQUENCE,STRAND,Error_Profile,Quality_Profile):
    new_sequence=''
    quality_score=''
    for i in range(len(OLD_SEQUENCE)):
        position=i+1
        new_sub,mutant_index=get_base(OLD_SEQUENCE[i],position,Error_Profile,STRAND)
        new_sequence=new_sequence+new_sub
        quality_score=quality_score+get_base_quality(new_sub,position,STRAND,mutant_index,Quality_Profile)
    return [new_sequence,quality_score]
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
def base_quality(Quality):
    Quality=Quality+round(random.normalvariate(0,2)) #a normal distribution with mean=0,sd=2
    if Quality<33:
        return chr(33)
    elif Quality>74:
        return chr(74)
    else:
        return chr(Quality)
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