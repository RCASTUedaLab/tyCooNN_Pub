import os
import inference.InferenceAndMakeBam as inference

def testEvaluate():


    flist = []

    indirs = ["/mnt/share/ueda/fast5/basecalled/200713_ecoli_TruB_total_v7/workspace"]
    outpath = "/mnt/share/ueda/TyCooNNPubTest/trub"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    flist.append((outpath,indirs))
    #
    indirs = ["/data/suzukilab/seqdata/basecall/210415_ecoli_ThiI_total_LB_BW_sta_V7/workspace",
              "/data/suzukilab/seqdata/basecall/210514_ecoli_ThiI_total_LB_sta_1_V7/workspace",
              "/data/suzukilab/seqdata/basecall/210519_ecoli_ThiI_total_LB_BW_sta_2_V7/workspace",
              "/data/suzukilab/seqdata/basecall/210519_ecoli_ThiI_total_LB_BW_sta_3_V7/workspace"]
    outpath = "/mnt/share/ueda/TyCooNNPubTest/thii"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    flist.append((outpath,indirs))




    indirs = ["/data/suzukilab/seqdata/basecall/201028_ecoli_WT_total_LB_BW_sta_2_V7/workspace",
              "/data/suzukilab/seqdata/basecall/210213_ecoli_WT_total_LB_BW_sta_1_V7/workspace",
              "/data/suzukilab/seqdata/basecall/210213_ecoli_WT_total_LB_BW_sta_2_V7/workspace"]
    outpath = "/mnt/share/ueda/TyCooNNPubTest/wt"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    flist.append((outpath,indirs))

    # indirs = "/mnt/share/ueda/fast5/basecalled/210514_ecoli_ThiI_total_LB_sta_1_V7/workspace/no_sample/20210514_0548_MN29778_FAP33126_bf303fcb/fast5"
    # outpath = "/mnt/share/ueda/TyCooNNPub/tapT"
    # if not os.path.exists(outpath):
    #     os.makedirs(outpath)
    # flist.append((outpath,indirs))

    #
    #
    # indirs = "/mnt/share/ueda/fast5/basecalled/201215_ecoli_TrmB_total_LB_BW_sta_V7/workspace"
    # outpath = "/mnt/share/ueda/TyCooNNPub/trmb"
    # if not os.path.exists(outpath):
    #     os.makedirs(outpath)
    # flist.append((outpath,indirs))




    # indirs = "/data/suzukilab/seqdata/basecall/201014_ecoli_WT_total_LB_BW_log_V7/workspace"
    # outpath = "/mnt/share/ueda/TyCooNNPub/log"
    # if not os.path.exists(outpath):
    #     os.makedirs(outpath)
    # flist.append((outpath,indirs))
    #
    # indirs = "/data/suzukilab/seqdata/basecall/201028_ecoli_WT_total_LB_BW_sta_1_V7/workspace"
    # outpath = "/mnt/share/ueda/TyCooNNPub/sta1"
    # if not os.path.exists(outpath):
    #     os.makedirs(outpath)
    # flist.append((outpath,indirs))


    paramPath = '/mnt/share/ueda/TyCooNNPub/setting.yaml'


    ref = "/mnt/share/ueda/trna_data/ref/ecolitRNA_unmod_full.fa"
    refdir = "/mnt/share/ueda/trna_data/ref/"
    # configdir = "/home/ueda/project/tyCooNN_Pub/resource/model_comb_ivtmap_6k_pub"
    configdir = "/mnt/share/bhaskar/combined_train_rcc_ivt/model_comb_ivtmap_rcc_10k_pub"

    for outpath, indirs in flist:

        # evaluate(paramPath, indirs, outdir, outpath, fasta, fasta5out, threshold=0.75)
        inference.evaluate(paramPath,indirs,configdir,outpath,ref,refdir,threshold=0.65)


import tensorflow as tf
with tf.device('/CPU:0'):
    testEvaluate()




