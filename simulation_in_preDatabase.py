import sys
sys.path.append(r"/data/users/liuyunhe/gitup")
import SSCRNA as sc
import MyProcess10
import time
import os
import pickle
from multiprocessing import Manager,Process,Pool
#
seq_size = int(sys.argv[1])
WorkPath = sys.argv[2]
Process = int(sys.argv[3])
Filepath=WorkPath+"Simulation_Path/"
os.chdir(Filepath)
#
with open(WorkPath+'cell_barcord_gene_name_seq.pickle', 'rb') as handle:
    xc=pickle.load(handle)
cell_barcord,f_na,f_seq=xc[0],xc[1],xc[2]
#
frag_list_path=(WorkPath+r'fragment_data/')
pcr_list_path=(WorkPath+r'pcr_control_data/')
umi_path=(WorkPath+r'umi_data/')
#
real_size,PCR_list=sc.get_search_key(pcr_list_path)
#
error_profile_r=sc.Read_Error(WorkPath+'Error_Profile.txt')
quality_profile=sc.Read_Quality(WorkPath+'Quality_Profile.txt')
#
lib_index=sc.get_library_index(real_size,seq_size)
count_sparse=sc.get_all_Count_UMI(pcr_list_path,PCR_list,lib_index)
lib_index_fragment=sc.get_fragment_seq2(frag_list_path,umi_path,count_sparse)
#leave it
adaptor=''
#
def get_chunck_mission(library_size,chunk_size):
    mission_pool=[]
    start_q=0
    end_q=chunk_size
    while end_q>start_q:
        mission_pool.append(range(start_q,end_q))
        start_q=start_q+chunk_size
        end_q=min(library_size,(end_q+chunk_size))
    return mission_pool

def list_to_str(list_a):
    mer_a=[]
    mer_a.append('[')
    for sub_a in list_a:
        mer_a.append(str(sub_a))
        mer_a.append(',')
    mer_a.pop()
    mer_a.append(']')
    return ''.join(mer_a)
#1
MP=get_chunck_mission(len(lib_index_fragment),int(seq_size/Process))
M_index=[]
i_index=0
for each in MP:
    M_index.append(i_index)
    i_index=i_index+1
###
while len(M_index)>0:
    #make mission queue
    workqueue=Manager().Queue()
    for i in MP[M_index[0]]:
        workqueue.put(i)
    size_queue=workqueue.qsize()
    #make process list
    process_list=[]
    for j in range(Process):
        process_list.append(MyProcess10.myThread(('Process_'+str(j)+'-Chunk_'+str(M_index[0])),workqueue,Filepath,'simulation_sequence',f_seq,f_na,adaptor,cell_barcord,error_profile_r,quality_profile,lib_index_fragment))
    #start
    try:
        for P in process_list:
            P.start()
        while not workqueue.empty():
            time.sleep(1)
            now_size=workqueue.qsize()
            print('\r'+'▇'*((((size_queue-now_size)*100)//(size_queue))//2)+str((((size_queue-now_size)*100)//(size_queue)))+'%', end='')
    except ConnectionResetError:
        print('\n'+'chunk '+str(M_index[0])+' raise error!!!')
        for P in process_list:
            P.terminate()
        os.system('rm '+Filepath+'simulation_sequence'+'_line_Process_*-Chunk_'+str(M_index[0])+'_R[1,2].fq')
    else:
        print('\n'+'chunk '+str(M_index[0])+'generate is complete!!!')
        for P in process_list:
            P.terminate()
        os.system('cat '+Filepath+'simulation_sequence'+'_line_Process_*-Chunk_'+str(M_index[0])+'_R1.fq >'+Filepath+'simulation_sequence'+'_line_'+str(M_index[0])+'_R1.fq')
        os.system('cat '+Filepath+'simulation_sequence'+'_line_Process_*-Chunk_'+str(M_index[0])+'_R2.fq >'+Filepath+'simulation_sequence'+'_line_'+str(M_index[0])+'_R2.fq')
        os.system('rm '+Filepath+'simulation_sequence'+'_line_Process_*-Chunk_'+str(M_index[0])+'_R[1,2].fq')
        print('\n'+'chunk '+str(M_index[0])+'merge is complete!!!')
        del(M_index[0])
#
for P in process_list:
    P.terminate()
