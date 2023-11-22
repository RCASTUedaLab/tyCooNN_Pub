import time
import numpy as np
from ont_fast5_api.fast5_interface import get_fast5_file

def get_fast5_reads_from_file(fast5_filepath:str):
    reads = []
    start_time = time.time()
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read in f5.get_reads():
            readid = read.read_id
            read_info = read.handle[read.raw_dataset_group_name].attrs
            mux = read_info['start_mux']
            duration = read_info['duration']
            median_before = read_info['median_before']
            read_start_time = read_info['start_time']
            channel_info = read.get_channel_info()
            channel_number = channel_info['channel_number']

            row = read.handle["Raw"]
            signal = row["Signal"][()]
            digitisation = channel_info['digitisation']
            offset = channel_info['offset']
            range_value = channel_info['range']
            pA_signal = (signal+offset)*range_value/digitisation

            analyses = read.list_analyses()
            basecall_run = None
            mapinfo_run = None
            for analysis,bpath in analyses:
                trace = read.get_analysis_dataset(bpath,"BaseCalled_template/Trace")
                if trace is not None: basecall_run = bpath
                if 'filterflg' in read.get_analysis_attributes(bpath): mapinfo_run = bpath
            if basecall_run is None:
                basecall_run = read.get_latest_analysis("Basecall_1D")
            fastq = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Fastq")
            trace = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Trace")
            move = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Move")

            if mapinfo_run is None:
                filterflg = -1
                filterpass = False
                trimSuccess = False
                tRNA_infer = None
                tRNAIndex = -1
                softmax_prob = 0.0
            else:
                attrs = read.get_analysis_attributes(mapinfo_run)
                filterflg = attrs['filterflg']
                filterpass = attrs['filterpass']
                trimSuccess = attrs['trimSuccess']
                tRNA_infer = attrs['tRNA']
                tRNAIndex = attrs['tRNAIndex']
                softmax_prob = attrs['value']
            map_attrs = {'filterflg': filterflg, 'filterpass': filterpass, 'trimSuccess': trimSuccess, 
                         'tRNA_infer': tRNA_infer, 'tRNAIndex': tRNAIndex, 'softmax_prob': softmax_prob}

            if len(trace) >0:
                read = Read(read_id=readid,signal=pA_signal,trace=trace,move=move,fastq=fastq,\
                    channel_number=channel_number,mux=mux,duration=duration,median_before=median_before,\
                    start_time=read_start_time,map_attrs=map_attrs)
                reads.append(read)
        print('       {}: {}reads {}s'.format(fast5_filepath,len(reads),time.time()-start_time)) # logging
        return reads


ascii_order = '!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
ascii_score_dict = {ascii_order[k]:k for k in range(len(ascii_order))}


class Read():
    """
    This class contains the all information about each read.
    The result of mapping and viterbi is added by using add_mapping_result() and add_viterbi_result()

    By using get_move_positions, you can get the segmentation data of guppy, which contains indexes of segmentation corresponding to signal.
    """
    def __init__(self,read_id,signal,trace,move,fastq,channel_number,mux,duration,median_before,start_time,map_attrs):
        self.read_id = read_id

        self.adapter_signal =  signal[:len(signal)-10*len(trace)][::-1].astype(np.float64)
        self.signal = signal[len(signal)-10*len(trace):][::-1].astype(np.float64)
        self.trace = trace[::-1].astype(np.int16)
        self.move = move[::-1].astype(np.int16)
        self.fastq = fastq

        fastq_list = self.fastq.split('\n')
        self.sequence = fastq_list[1].replace('T','U')
        self.qscore = np.array([ascii_score_dict[symbol] for symbol in fastq_list[3]],dtype=np.int16)
        self.mean_qscore = sum(self.qscore)/len(self.qscore)
        self.channel_number = channel_number
        self.mux = mux
        self.duration = duration
        self.median_before = median_before
        self.start_time = start_time

        self.map_attrs = map_attrs

        self.normalized_signal = None

        self.mapping_results = []
        self.viterbi_results = []
        self.processid = 0

        #used for trim and normalize
        self.trimIdxbyHMM = 0
        self.trimIdxbyMapping = 0
        self.normalizeDelta = 0
        self.trimmedSignal = []
        self.trimmedTrace = []
        self.countcpd = 0
        self.normalizemed = 0

        self.formatSignal = []

        #trim
        self.trimSuccess = False
        self.filterFlg = 0

    def normalize(self,shift,scale):
        self.normalized_signal = (self.signal - shift) / scale
        
    def setProcessId(self,n):
        self.processid = n
    def getProcessId(self):
        return self.processid

    def add_summary_information(self,channel,mux,start_time,duration,template_start,template_duration):
        self.channel = channel
        self.mux = mux
        self.start_time = start_time
        self.duration = duration
        self.template_start = template_start
        self.template_duration = template_duration
    
    def get_move_positions(self):
        """
        return the list of segmentation positions by guppy. [(index,base),...]
        """
        positions = []
        seq_num = 0
        for n,m in enumerate(self.move):
            index = (n+1)*10
            if m ==1:
                seq_num += 1
                if seq_num >= len(self.sequence):
                    break
                positions.append((index,self.sequence[seq_num]))
        return positions
    
    def __repr__(self):
        print_str = ''
        print_str += '\nRead id:{}\n'.format(self.read_id)
        print_str += ' Signal Length:{}\n'.format(len(self.signal))
        print_str += ' Sequence:{}\n'.format(self.sequence)
        print_str += '   qscore:{}\n'.format(''.join([str(score) + ' ' for score in self.qscore]))
        print_str += ' Sequence Length:{}\n'.format(len(self.sequence))
        print_str += ' Channel:{} Mux:{}\n'.format(self.channel,self.mux)
        print_str += ' Whole duration:{}, Template duration:{}\n'.format(self.duration,self.template_duration)
        return print_str
