import sys
sys.path.append("../")
import yaml
import inference.Inference as inference

import os

def testEvaluate(opts):
    inference.evaluate(opts)

input_options = "inp_opt.yaml"
dbname = sys.argv[1]
given_name = sys.argv[2]

io = {'i':"/data/suzukilab/seqdata/basecall/",
      'o':"/mnt/share/bhaskar/rcc_data_extend/infer_98/"}

with open(input_options) as f:
    opts = yaml.safe_load(f)
    opts['inp_loc'] = io['i'] + dbname + "/workspace"
    opts['out_loc'] = io['o'] + given_name
    testEvaluate(opts) 
