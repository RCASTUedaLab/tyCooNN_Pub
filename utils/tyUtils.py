import glob
from utils.tyRead import Read
from utils.tyParam import tyParam
import multiprocessing
from multiprocessing import Pool
from ont_fast5_api.fast5_interface import get_fast5_file
import yaml, sys, itertools
import matplotlib.pyplot as plt
import csv

def get_parameter(yaml_path:str):
    print('loading parameter from {}\n'.format(yaml_path))
    try:
        with open(yaml_path) as file:
            return tyParam(yaml.safe_load(file))

    except Exception as e:
        print('Exception occurred while loading YAML...', file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)

def split_list(l, n=1000):

    for idx in range(0, len(l), n):
        yield l[idx:idx + n]

def get_fast5_reads(directory:str,MAX_CORE:int,readmax = -1):
    """
    return the list of reads from fast5 files in the directory

    Args:
        directory (str): path to the directory containing fast5 files
        MAX_CORE (int): maximum number of cores

    Returns:
        list: list of Read intances

    """
    print('load fast5 reads from {}'.format(directory))
    f5list =  get_fast5_files_in_dir(directory,MAX_CORE)
    #print('fast5 list {}'.format(f5list))
    if len(f5list) == 1:
       reads = get_fast5_reads_from_file(f5list[0])
       print('Finish. 1 reads are loaded\n'.format(len(reads)))
       return reads
    if readmax > 0:
        upto = min(readmax,len(f5list)-1)
        f5list = f5list[0:upto]
    ncore = get_number_of_core(MAX_CORE=MAX_CORE)
    with Pool(ncore) as p:
        reads = p.map(get_fast5_reads_from_file, f5list)
        reads = list(itertools.chain.from_iterable(reads)) # flatteining
    print('Finish. {}reads are loaded\n'.format(len(reads)))
    return reads

def get_fast5_reads_dirs(directories:list,MAX_CORE:int,readmax = -1):
    """
    return the list of reads from fast5 files in the directory

    Args:
        directories (str): path to the directory containing fast5 files
        MAX_CORE (int): maximum number of cores

    Returns:
        list: list of Read intances

    """
    f5list = []
    for directory in directories:
        print('load fast5 reads from {}'.format(directory))
        f5list.extend(get_fast5_files_in_dir(directory,MAX_CORE))

    #print('fast5 list {}'.format(f5list))
    if len(f5list) == 1:
       reads = get_fast5_reads_from_file(f5list[0])
       print('Finish. 1 reads are loaded\n'.format(len(reads)))
       return reads
    if readmax > 0:
        upto = min(readmax,len(f5list)-1)
        f5list = f5list[0:upto]
    ncore = get_number_of_core(MAX_CORE=MAX_CORE)
    with Pool(ncore) as p:
        reads = p.map(get_fast5_reads_from_file, f5list)
        reads = list(itertools.chain.from_iterable(reads)) # flatteining
    print('Finish. {}reads are loaded\n'.format(len(reads)))
    return reads

def get_number_of_core(MAX_CORE:int):
    ncore = multiprocessing.cpu_count()
    if ncore > MAX_CORE:
        ncore = MAX_CORE
    return ncore


def get_fast5_files_in_dir(directory:str,MAX_CORE:int):
    ncore = get_number_of_core(MAX_CORE=MAX_CORE)
    filelist = list(sorted(glob.glob(directory + '/**/*.fast5',recursive=True)))
    print("Number of fast5 files: %d" % len(filelist))
    with Pool(ncore) as p:
        ok_test = p.map(check_ok_file,filelist)
    filelist2 = []
    for i,file in enumerate(filelist):
        if ok_test[i]:
            filelist2.append(file)
    print("Number of fast5 files [ok]: %d" % len(filelist2))
    return filelist2

def check_ok_file(file:str):
    # print("check %s" % file)
    try:
        with get_fast5_file(file,mode='r') as f5:
            for read in f5.get_reads():
                readid = read.read_id
                row = read.handle["Raw"]
                signal = row["Signal"][()]
                channel_info = read.get_channel_info()
                digitisation = channel_info['digitisation']
                offset = channel_info['offset']
                range_value = channel_info['range']
                read_info = read.handle[read.raw_dataset_group_name].attrs
                duration = read_info['duration']
                basecall_run = read.get_latest_analysis("Basecall_1D")
                fastq = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Fastq")
                trace = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Trace")
                move = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Move")
            ok = True
    except:

        ok = False
    print(file,ok)
    return ok

def getOrNone(groups,partialKey):
    r = None
    if groups is None:
        return None
    for k in groups:
        if partialKey in k:
            r = groups[k]
    return r

def get_fast5_reads_from_file(fast5_filepath:str):
    reads = []
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read in f5.get_reads():
            readid = read.read_id

            row = read.handle["Raw"]
            signal = row["Signal"][()]
            channel_info = read.get_channel_info()
            digitisation = channel_info['digitisation']
            offset = channel_info['offset']
            range_value = channel_info['range']
            pA_signal = (signal + offset) * range_value / digitisation
            read_info = read.handle[read.raw_dataset_group_name].attrs
            duration = read_info['duration']

            # Search for latest analysis dataset
            basecall_run = read.get_latest_analysis("Basecall_1D")

            fastq = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Fastq")
            trace = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Trace")
            move = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Move")

            #read_id, signal, tracelen, fastq
            if len(trace) >0:
                #read_id, signal, trace, move, fastq, duration):
                read = Read(read_id=readid,signal=pA_signal,trace=trace,move=move,fastq=fastq,duration=duration)
                reads.append(read)

    return reads


def plot_history(history,
                 save_graph_img_path,
                 fig_size_width,
                 fig_size_height,
                 lim_font_size):
    acc = history.history['accuracy']
    val_acc = history.history['val_accuracy']
    loss = history.history['loss']
    val_loss = history.history['val_loss']

    epochs = range(len(acc))

    # Ot\
    plt.figure(figsize=(fig_size_width, fig_size_height))
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = lim_font_size  # StHg
    # plt.subplot(121)

    # plot accuracy values
    plt.plot(epochs, acc, color="blue", linestyle="solid", label='train acc')
    plt.plot(epochs, val_acc, color="green", linestyle="solid", label='valid acc')
    # plt.title('Training and Validation acc')
    # plt.grid()
    # plt.legend()

    # plot loss values
    # plt.subplot(122)
    plt.plot(epochs, loss, color="red", linestyle="solid", label='train loss')
    plt.plot(epochs, val_loss, color="orange", linestyle="solid", label='valid loss')
    # plt.title('Training and Validation loss')
    plt.legend()
    plt.grid()

    plt.savefig(save_graph_img_path)
    plt.close()  # obt@

def printResult(out, ep, model, test_x):
    prediction = model.predict(test_x)
    outFile = out + "/result" + str(ep) + '.csv'
    with open(outFile,'w') as fw:
        writer = csv.writer(fw, lineterminator='\n')
        for row in prediction:
            writer.writerow(row)
