import multiprocessing
from S_S2 import get_full_seq2,get_one_read,get_seq_read,get_seq_title2
class myThread (multiprocessing.Process):
    def __init__(self, name, q,FILE_PATH,FILE_NAME,SEQUENCES,SEQUENCE_NAMES,ADAPTOR,Cell_Barcord,ERROR_Profile,QUALITY_Profile,LIBRARY_INDEX_Iterm):
        super().__init__()
        self.name=name
        self.q = q
        self.FILE_PATH = FILE_PATH
        self.FILE_NAME = FILE_NAME
        self.SEQUENCES = SEQUENCES
        self.SEQUENCE_NAMES = SEQUENCE_NAMES
        self.ADAPTOR = ADAPTOR
        self.Cell_Barcord = Cell_Barcord
        self.ERROR_Profile = ERROR_Profile
        self.QUALITY_Profile = QUALITY_Profile
        self.LIBRARY_INDEX_Iterm = LIBRARY_INDEX_Iterm
    def run(self):
        print (self.name+' start working now!')
        process_data(self.q,self.FILE_PATH,self.FILE_NAME,self.SEQUENCES,self.SEQUENCE_NAMES,self.ADAPTOR,self.Cell_Barcord,self.ERROR_Profile,self.QUALITY_Profile,self.LIBRARY_INDEX_Iterm) #in quene, will not stop
        #print (self.name+' have finish his work!')
#add 16.2
def process_data(q,FILE_PATH,FILE_NAME,SEQUENCES,SEQUENCE_NAMES,ADAPTOR,Cell_Barcord,ERROR_Profile,QUALITY_Profile,LIBRARY_INDEX_Iterm): 
    FILE=FILE_PATH+'/'+FILE_NAME+'_line_'+multiprocessing.current_process().name
    flush_empty_index=0
    write_flush_1=[]
    write_flush_2=[]
    while 1:
        if not q.empty():
            fragment_bed_full = LIBRARY_INDEX_Iterm[q.get()]
            full_seq=get_full_seq2(SEQUENCES,ADAPTOR,Cell_Barcord,fragment_bed_full)
            read1_seq=get_one_read(full_seq,1,[0,1],125)
            read1_all=get_seq_read(read1_seq,1,ERROR_Profile,QUALITY_Profile)
            read2_seq=get_one_read(full_seq,2,[0,1],125)
            read2_all=get_seq_read(read2_seq,2,ERROR_Profile,QUALITY_Profile)
            #insert to flush
            write_flush_1.extend([get_seq_title2(SEQUENCE_NAMES,Cell_Barcord,fragment_bed_full,1),'\n'])
            write_flush_1.extend([read1_all[0],'\n'])
            write_flush_1.extend(['+','\n'])
            write_flush_1.extend([read1_all[1],'\n'])
            write_flush_2.extend([get_seq_title2(SEQUENCE_NAMES,Cell_Barcord,fragment_bed_full,2),'\n'])
            write_flush_2.extend([read2_all[0],'\n'])
            write_flush_2.extend(['+','\n'])
            write_flush_2.extend([read2_all[1],'\n'])
            flush_empty_index=1
            if len(write_flush_1)<1000:
                pass
            else:
                with open(FILE+'_R1.fq','a') as f1,open(FILE+'_R2.fq','a') as f2:
                    f1.writelines(write_flush_1)
                    f2.writelines(write_flush_2)
                write_flush_1=[]
                write_flush_2=[]
                flush_empty_index=0
        elif flush_empty_index==1:
            with open(FILE+'_R1.fq','a') as f1,open(FILE+'_R2.fq','a') as f2:
                f1.writelines(write_flush_1)
                f2.writelines(write_flush_2)
            write_flush_1=[]
            write_flush_2=[]
            flush_empty_index=0