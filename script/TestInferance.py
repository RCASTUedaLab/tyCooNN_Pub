import sys
sys.path.append("../")
import yaml
import inference.Inference as inference

import os

def testEvaluate(opts):
    '''
    opt entries:
    inp_loc = "/share/trna/testdata/201026_ecoli_WT_mcmo5U_total_LB_BW_sta_V7_workspace/"
    out_loc = "/home/bhaskar/work/tyCooNN_2202_master/data/out/WT_mcmo5U_total_LB_BW_sta"

    model_loc = "/home/bhaskar/work/tyCooNN_2202_master/data/out"
    fasta_loc = "/share/bhaskar/tyCooNNTest/fasta/ecolitRNA_full.fa"
    param_loc = "/home/bhaskar/work/tyCooNN_2202_master/setting.yaml"
    post_filter_threshold = 0.75
    '''

    inference.evaluate(opts)

input_options = sys.argv[1]
with open(input_options) as f:
    opts = yaml.safe_load(f)
    testEvaluate(opts) 
