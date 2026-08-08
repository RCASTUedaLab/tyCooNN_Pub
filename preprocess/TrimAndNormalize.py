from multiprocessing.pool import Pool

import statistics
import numpy as np
from scipy import ndimage as ndi
import copy
from hmmlearn import hmm
from Bio import pairwise2
from numba import jit
from functools import partial
import random

def _down_sampling(a):
    if (len(a) % 2 == 1): a = np.append(a, 0)
    c = copy.copy(a)
    c = np.array(c)
    ret = c.reshape(-1, 2).mean(axis=1)

    return ret

def down_sampling(a, deg):
    for i in range(deg):
        a = _down_sampling(a)

    return a


def getStartIndexes(signal, downSampleDegree=4, deltaThreshold=20, lenThreshold=500, minSiglen=2000,maxSiglen=10000):

    dssig = down_sampling(signal, downSampleDegree)
    diff = np.array([1, 0, -1])
    diffsig = ndi.convolve(dssig, diff)
    fsig = np.where(diffsig > deltaThreshold)[0]

    factor = (2 ** 4)
    dd = fsig * factor
    data = []
    minO = minSiglen
    maxO = maxSiglen
    last = 0
    for i in range(len(dd)):
        v = dd[i]
        if (v > minO and v < maxO and (v - last) > lenThreshold):
            data.append(v)
        last = v
    return data

def applyHMM(signal,minIdx=500):

    n_components = 2
    model = hmm.GaussianHMM(n_components=n_components)
    model.startprob_ = np.array([1.0, 0.0])

    model.transmat_ = np.array([[0.5, 0.5],
                                [0, 1]])

    model.means_ = np.array([[0.66],
                             [-0.45]])

    model.covars_ = np.array([[0.2],
                              [0.1]])

    signallen = len(signal)
    #sub = signal[0:signallen - 1100]
    sub = signal
    sub = sub.reshape(-1, 1)
    #print(sub)
    r = model.predict(sub)
    #print(r)
    end = 0

    if (np.any(r == 1)):
       end = np.where(r == 1)[0][0]
    #     atleast 2100 signal point    after adaptor
    #     not longer than 5600
    leflen = signallen - end
    if (end < minIdx):
        return 0

    return end


@jit(nopython=True)
def get_start_and_end_index(query_sequence,ref_sequence,query_length,ref_length):

    q_idx = -1
    r_idx = -1
    r_st = -1
    r_en = -1
    q_st = -1
    q_en = -1
    for n in range(len(query_sequence)):
        if query_sequence[n] != '-':
            q_idx += 1
        if ref_sequence[n] != '-':
            r_idx += 1
            if r_st == -1 and query_sequence[n] == ref_sequence[n]:
                # alignment start
                r_st = r_idx
                q_st = q_idx
        if q_idx+1 == query_length or r_idx+1 == ref_length:
            r_en = r_idx
            q_en = q_idx
            break
    return np.array([r_st,r_en,q_st,q_en],dtype=np.int64)

def applySeq(ref, read):

    fasta = read.sequence[::-1] # to from 3P
    ref = ref.replace('U','T')[::-1] #
    hang = 5
    fasta = fasta[0:len(ref) + hang]
    match,mismatch,indel,extension = 16, -12, -30, -15
    alignment = pairwise2.align.localms(fasta, ref, match,mismatch,indel,extension)
    query_sequence = alignment[0].seqA
    ref_sequence = alignment[0].seqB
    r_st, r_en, q_st, q_en = get_start_and_end_index(query_sequence, ref_sequence, len(fasta),
                                                     len(ref))
    moveform3p = read.move[::-1] # to from 3P
    return getBound(moveform3p, q_en)

def getBound(moveform3p, q_en):

    cnt = 0
    nucidx = 0
    for flg in moveform3p:

        if flg > 0:
            nucidx += 1
        if nucidx == q_en:
            return cnt * 10
        cnt += 1
    return 0

def getHighAGPeak(sigfrom3p):

    m = 0
    medmax = 0
    for n in range(5):
        firstsigfrom3p = sigfrom3p[m:m + 500]
        if len(firstsigfrom3p) == 0:
            break
        med = statistics.median(firstsigfrom3p)
        if medmax < med:
            medmax = med
        m = m + 100
    return medmax


def filterFlg(read,param,trimSuccess):

    if read.mean_qscore < param.qval_min:
        return 1

    siglen = len(read.signal)
    if siglen > param.signallen_max:
        return 2

    if (read.duration / siglen) < param.duratio_rate_max:
        return 3

    if read.normalizeDelta < param.delta_min:
        return 4

    if read.normalizeDelta > param.delta_max:
        return 5

    readlen = len(read.sequence)
    if readlen <= param.readlen_min:
        return 6

    if readlen >= param.readlen_max:
        return 7

    if not trimSuccess:
        return 8

    return 0

def trimAdaptorEach(read,param):

    sigfrom3p = read.signal
    adaptor1med = 0
    if len(read.adapter_signal) > 0:
        adaptor1med = statistics.median(read.adapter_signal)

    # take median of first 50-550
    adaptor2med = getHighAGPeak(sigfrom3p)
    diffadop = adaptor2med - adaptor1med
    hmmIndex = 0
    normsig = None
    mediantoset = 0
    if diffadop > 20:
        diffthery = param.adap2thery - param.adap1thery
        step = diffadop / diffthery
        mediantoset = adaptor1med + (step * (param.meantoSet - param.adap1thery))
        normsig = (sigfrom3p - mediantoset) / diffadop
        hmmIndex = applyHMM(normsig)


    mappingIndex = applySeq(param.firstAdaptor, read)

    read.trimIdxbyHMM = hmmIndex
    read.trimIdxbyMapping = mappingIndex
    read.normalizeDelta = diffadop
    trimSafeIdx = min(hmmIndex, mappingIndex)
    if trimSafeIdx == 0:
        trimSafeIdx = max(hmmIndex, mappingIndex)
    #
    trimNormSig = None
    if normsig is not None:
        trimNormSig = normsig[trimSafeIdx:len(normsig)]

    read.trimmedSignal = trimNormSig
    read.normalizemed = mediantoset

    read.trimSuccess = trimSafeIdx > 0 and trimNormSig is not None
    read.filterFlg = filterFlg(read,param,read.trimSuccess)

    return read

#trim
def trimAdaptor(reads,param):

    trimFunction = partial(trimAdaptorEach, param=param)
    with Pool(param.ncore) as p:
        reads = p.map(trimFunction, reads)
    return reads

#fromat to length
def formatSignal(reads, param):

    formatFunction = partial(_format, param=param)
    with Pool(param.ncore) as p:
         reads = p.map(formatFunction, reads)
    return reads

def _format(read,param):

    trimsignal = read.trimmedSignal
    if read.trimSuccess:
        formatsig = binned(trimsignal, param.trimlen)
        zoformatsig = zeroToOne(formatsig)
        read.formatSignal = zoformatsig
    return read

# Edit Nov 8, 2022
def zeroToOne(formatsig):


    formatsig = np.clip(formatsig, -1.5, 1.5)
    formatsig = formatsig + 1.5  / 3
    return np.clip(formatsig, 0, 1)

#def zeroToOne(formatsig):
#
#    a = formatsig + 1 / 2
#    return np.clip(a, 0, 1)

def binned(trimsignal, trimlength, mode=1):

    if len(trimsignal) == trimlength:
        return trimsignal  # not very likely to happen

    if len(trimsignal) > trimlength:
        # trim from first
        return trimsignal[0:trimlength]
    else:
        #
        ret = np.zeros(trimlength)
        diff = np.array([1, 0, -1])

        med = statistics.median(trimsignal)
        diffsig = ndi.convolve(trimsignal, diff)
        sigma = np.std(diffsig) / 10

        siglen = len(trimsignal)
        left = trimlength - siglen
        lefthalf = left // 2
        #
        for n in range(trimlength):
            if n < lefthalf or n >= (trimlength - lefthalf):
                if mode == 1:
                    ret[n] = med + noise(sigma)
                else:
                    ret[n] = 0  # zero pad
            else:
                idxa = n - lefthalf - 1
                if idxa < 0:
                    idxa = 0
                ret[n] = trimsignal[idxa]
        return ret

@jit(nopython=True)
def noise(sigma):

    return random.gauss(0, sigma)





